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

import java.util.List;
import java.util.Collections;
import java.util.Iterator;
import java.util.ListIterator;

/**
 * @author Kerstin Scheubert
 */
public class Spectrum  implements java.io.Serializable, Comparable<Spectrum> {

	private static final long serialVersionUID = 4085155682101803641L;
	private int collisionEnergy;
	private List<Peak> peaks;
	private double tic;
	private double maxIntensity;
	private double minIntensity;
	private double RT;
	private int mslevel;
	private double parentMass;
	private Peak parentPeak;

	public Spectrum(int collisionEnergy, double tic, List<Peak> peaks){
		this (collisionEnergy, tic, peaks, 0.0);
	}

	public Spectrum(int collisionEnergy, List<Peak> peaks){
		this(collisionEnergy, 0.0, peaks);
	}	

	public Spectrum(int collisionEnergy, double tic, List<Peak> peaks, double RT) {
		this.collisionEnergy = collisionEnergy;
		this.tic = tic;
		this.peaks = peaks;
		this.RT = RT;
		findMinMaxInt();    
		this.mslevel=2;
	}
	
	public Spectrum(int collisionEnergy, double tic, List<Peak> peaks, double parentmass, int mslevel){
		this(collisionEnergy, tic, peaks, 0.0, parentmass, mslevel);
	}

	public Spectrum(int collisionEnergy, List<Peak> peaks, double parentmass, int mslevel){
		this(collisionEnergy, 0.0, peaks, 0.0, parentmass, mslevel);
	}	

	public Spectrum(int collisionEnergy, double tic, List<Peak> peaks, double RT, double parentmass, int mslevel) {
		this(collisionEnergy, tic, peaks, RT);
		for (Peak p : peaks) {			
			p.addSpectrum(this);			
		}
		this.mslevel=mslevel;
		this.parentMass=parentmass;		
	}

	public int getCollisionEnergy(){
		return collisionEnergy;
	}

	public double getTic(){
		return tic;
	}

	public List<Peak> getPeaks(){
		return peaks;
	}

	public void setCollisionEnergy(int collisionEnergy){
		this.collisionEnergy = collisionEnergy;
	}

	public void setTic(double tic){
		this.tic = tic;
	}

	public void setPeaks(List<Peak> peaks){
		this.peaks = peaks;
	}

	public double getRT() {
		return RT;
	}

	public void setRT(double rt) {
		RT = rt;
	}
	
	public int getMSLevel(){
		return mslevel;
	}
	
	public void setMSLevel(int mslevel){
		this.mslevel=mslevel;
	}
	
	public double getParentMass(){
		return parentMass;
	}
	
	public void setParentMass(double parentmass){
		this.parentMass=parentmass;
	}
	
	public Peak getParentPeak(){
		return parentPeak;
	}
	
	public void setParentPeak(Peak parentPeak){
		this.parentPeak=parentPeak;
	}

	public int compareTo(Spectrum another){
		return (collisionEnergy - another.collisionEnergy);
	}

	public void calculateSNR(){
		Collections.sort(peaks, Peak.intensityComparator());
		double number = peaks.size(), i = 1.0;
		//double prob = 1.0;

		for (Peak peak : peaks){
			//prob = lambda*Math.exp(-1*lambda*i/number);
			peak.setRelIntensity(i/number);
			++i;
		}
		findMinMaxInt();
	}

	public void strictFilter(double threshold){
		for (Iterator<Peak> pIter = peaks.iterator(); pIter.hasNext();){
			if (pIter.next().getRelIntensity() < threshold){
				pIter.remove();
			}
		}
		findMinMaxInt();
	}

	public void calculateRaw(){
		if (tic != 0.0) {
			double sum = 0.0;
			for (Peak peak : peaks){
				sum += peak.getIntensity();
			}		
			for (Peak peak : peaks){
				peak.setRelIntensity(tic*peak.getIntensity()/sum);
			}
		}
		findMinMaxInt();
	}

	public void calcRelIntensities(){
		double totalIntensity = 0.0;
		for (Peak peak : peaks){
			totalIntensity += peak.getIntensity();
		}
		for (Peak peak : peaks){
			peak.setRelIntensity(peak.getIntensity()*100.0/totalIntensity);
		}
		findMinMaxInt();
	}

	public void copyIntensities(){
		for (Peak peak : peaks){
			peak.setRelIntensity(peak.getIntensity());
		}
		findMinMaxInt();
	}

	public void findIsotopePeaks(Parameters param){

		Peak curr, mono = null, prev = null;	
		int isotopePeaks = 0;
		for (Iterator<Peak> pIter = peaks.iterator(); pIter.hasNext();) {
			curr = pIter.next();
			double diff = (prev == null)?Double.POSITIVE_INFINITY:curr.getMass() - prev.getMass();
			double error = param.getErrorForMass(curr.getMass())*10/*Math.sqrt(2)*/;
			//System.out.println("diff "+diff+" "+error);
			if (diff >= Element.NEUTRON_MASS-error && diff <= Element.NEUTRON_MASS+error) {
				if (mono != null && mono.getRelIntensity() > curr.getRelIntensity()){ // Only Isotope Pattern if intensity relations fit.
					mono.addIsotope(curr);
					//System.out.println("Peak "+curr.getMass()+" is an isotope");
					++isotopePeaks;
					pIter.remove();
				} else {
					curr.setOverlaps(true);
					mono = curr;
				}
			} else {
				mono = curr;	
			}
			prev = curr;
		}
		System.out.println("Found "+isotopePeaks+" isotope peaks");
	}


	public double getMinIntensity() {
		return minIntensity;
	}

	public double errorCorrFactor(double intensity, Parameters param) {
		/*double gradient = (1-param.errorIncrease)/(maxIntensity-minIntensity),
  				 offset = 1-gradient*maxIntensity;
  	return gradient*intensity+offset;*/
		return (-param.errorIncrease+1)*intensity/maxIntensity+param.errorIncrease;
	}

	public void findMinMaxInt() {
		maxIntensity = Double.NEGATIVE_INFINITY;
		minIntensity = Double.POSITIVE_INFINITY;
		for (Peak p : peaks){
			if (p.getRelIntensity() > maxIntensity){
				maxIntensity = p.getRelIntensity();
			}
			if (p.getRelIntensity() < minIntensity){
				minIntensity = p.getRelIntensity();
			}
		}

	}

	public void decomposePeaks(Decomposer decomp, Parameters param, int mode) {
		for (Peak p : peaks){
			p.decompose(decomp, param, mode);
		}
	}

	public void scoreDecompositions(Parameters param, int charge) {
		for (Peak p : peaks){
			p.scoreDecompositions(param, this, charge);
		}
	}

	public void filterDecompositions() {
		for (Peak p : peaks){
			p.filterDecompositions();
		}
	}

	public void mergePeaks(Parameters param){
		Collections.sort(peaks);
		ListIterator<Peak> left = peaks.listIterator(), right = null;
		while (left.hasNext()){
			Peak leftPeak = left.next(), mergePeak = leftPeak; //
			right = peaks.listIterator(left.nextIndex()); //right: one beyond the peaks to merge.
			while (right.hasNext() && Math.abs(mergePeak.getMass()-(mergePeak=right.next()).getMass()) < param.mergeThreshold){
			}
			int numberMerges = right.previousIndex()-left.previousIndex()-1; // one merge less than there are Peaks to be merged
			for(int i=0; i<numberMerges && left.hasNext(); ++i){
				mergePeak = left.next();
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
				leftPeak.mergeDecompositions(mergePeak, param);
				left.remove();
			}
		} 
	}

	public void normIntensities(double desiredMax) {
		double max = 0.0;
		for (Peak p : peaks){
			if (p.getIntensity() > max){
				max = p.getRelIntensity();
			}
		}
		double scaling = desiredMax/max;
		for (Peak p : peaks){
			p.setRelIntensity(p.getRelIntensity()*scaling);
		}

		findMinMaxInt();
	}

}
