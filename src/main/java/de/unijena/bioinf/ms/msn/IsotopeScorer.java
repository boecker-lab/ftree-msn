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

import java.util.Iterator;
import java.util.List;

/**
 * @author Kerstin Scheubert
 */
public class IsotopeScorer {

  private Parameters param;
  private Spectrum spectrum;
  private List<Peak> sampleDistribution;
  private double scaling;
  
  public IsotopeScorer(Parameters param, Spectrum spectrum) {
    this.spectrum = spectrum;
    this.param = param;
  }
  
  public IsotopeScorer(List<Peak> sampleDistribution, Parameters param, Spectrum spectrum){
    this(param, spectrum);
    this.sampleDistribution = sampleDistribution;
    scaling = calculateScaling();
  }
  
  public List<Peak> getSampleDistribution() {
    return sampleDistribution;
  }

  public void setSampleDistribution(List<Peak> sampleDistribution) {
    this.sampleDistribution = sampleDistribution;
    scaling = calculateScaling();
  }

  public Parameters getParam() {
		return param;
	}

	public double getScaling() {
		return scaling;
	}

	public double getScore(List<Peak> candidateDistribution) {
    
    if (sampleDistribution.isEmpty() && candidateDistribution.isEmpty()) {
      return 0.0; // No changes to original score
    }
    
    double score = 0.0;
    final double sqrt2 = Math.sqrt(2);
    double diff = 0.0, sd = 0.0, error = 0.0;
    
    // scores masses of peaks
    Iterator<Peak> sampleIter = sampleDistribution.iterator(),
                   candidateIter = candidateDistribution.iterator();
    Peak sample = sampleIter.next(), candidate = candidateIter.next(); // monoisotopic peak is scored elsewhere
    Peak prevSample = null, prevCand= null;
    //Peak sample = null, candidate = null;
    
    while(sampleIter.hasNext() && candidateIter.hasNext()) {
      prevSample = sample;
      prevCand = candidate;
      sample = sampleIter.next();
      candidate = candidateIter.next();
      // score differences between the peaks => constant offset removed 
      diff = Math.abs(sample.getMass()-prevSample.getMass() - candidate.getMass() + prevCand.getMass());
      error = param.getErrorForMass(sample.getMass());
      sd = (error*spectrum.errorCorrFactor(sample.getRelIntensity(), param))/param.massDeviationPenalty;
      score = logAndAddScore(MathUtils.erfc(diff / (sd * sqrt2)), score);
    }
    
    // scores abundances of peaks
    sampleIter = sampleDistribution.iterator(); // reset Iterators
    candidateIter = candidateDistribution.iterator();
    while(sampleIter.hasNext() && candidateIter.hasNext()) {
      sample = sampleIter.next();
      candidate = candidateIter.next();
      diff = Math.abs(scaledIntensity(sample.getRelIntensity()) - candidate.getIntensity());
      error = param.intError*1e-2;//*sample.getRelIntensity();
      sd = error*spectrum.errorCorrFactor(sample.getRelIntensity(), param);
      score = logAndAddScore(MathUtils.erfc(diff / (sd * sqrt2)), score);
    }
    
    // score non-existence of last peak
    if (candidateIter.hasNext()){
      candidate = candidateIter.next();
      double sampleInt = spectrum.getMinIntensity();
      // scaling is determined for existing peaks only, thus the intensity of the imaginary peak is added to the scaling factor.  
      diff = candidate.getIntensity() - sampleInt/(scaling+sampleInt);
      if (diff > 0) {
        error = param.intError*1e-2;//*sampleInt;
        sd = error*spectrum.errorCorrFactor(sampleInt, param);
        score = logAndAddScore(MathUtils.erfc(diff / (sd * sqrt2))/2, score); //One sided erfc = erfc/2
      }
    }
    
    return score;
  }
  
  private static double logAndAddScore(double toAdd, double score) {
  	//return score + Math.log(toAdd);
		return score + 1 + Math.log(toAdd);
  	//return score * toAdd;
	}

	private double calculateScaling() {
    double intSum = 0.0;
    for (Peak p : sampleDistribution){
      intSum += p.getRelIntensity();
    }
    return intSum;
  }
  
	private double scaledIntensity(double intensity){
		// calculate scaled int, add offset and rescale by the second term
		return (intensity/scaling+param.intOffset)/(1+sampleDistribution.size()*param.intOffset);
	}
}