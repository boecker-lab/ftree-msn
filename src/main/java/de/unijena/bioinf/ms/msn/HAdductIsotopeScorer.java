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
import java.util.List;
import java.util.Iterator;

/**
 * @author Kerstin Scheubert
 */
public class HAdductIsotopeScorer extends IsotopeScorer{

  public HAdductIsotopeScorer(Parameters param, Spectrum spectrum) {
    super(param, spectrum);
  }
  
  public HAdductIsotopeScorer(List<Peak> sampleDistribution, Parameters param, Spectrum spectrum){
    super(sampleDistribution, param, spectrum);
  }
    
  public double getScore(List<Peak> candidateDistribution) {
    // Assertion    
    if (getSampleDistribution().size() != candidateDistribution.size()){
      System.out.println("Error size mismatch: "+getSampleDistribution().size()+" "+candidateDistribution.size());
    }
    
    List<Peak> shiftedDistribution = shift(candidateDistribution);
    double numerator = 0.0, denominator = 0.0;
    Iterator<Peak> sIter = shiftedDistribution.iterator(),
                   mIter = getSampleDistribution().iterator(); // measured peaks
    double sInt, mInt, csDiff; // Intensities, Difference
    for (Peak c : candidateDistribution){
      sInt = sIter.next().getIntensity();
      mInt = mIter.next().getRelIntensityDB(); // Normalized Intensities of measured spectra are in RelIntensityDB.
      csDiff = c.getIntensity()-sInt;
      numerator += (mInt-sInt)*csDiff;
      denominator += csDiff*csDiff;
    }
    double alpha = numerator / denominator;
    /*if (candidateDistribution.get(0).getMass() > 500 && candidateDistribution.get(0).getMass() < 550){
      System.out.println("alpha = "+alpha);
    }*/
    List<Peak> mixedDistribution;
    if (alpha > 1) {
      mixedDistribution = candidateDistribution;
    } else if (alpha < 0){
      mixedDistribution = shiftedDistribution;
    } else {
      mixedDistribution = mix(candidateDistribution, shiftedDistribution, alpha);
    }
    return super.getScore(mixedDistribution);
  }

  private List<Peak> shift(List<Peak> candidate){
    List<Peak> result = new ArrayList<Peak>(candidate.size()+1);
    result.add(new Peak(candidate.get(0).getMass(),0.0, 0));
    for (Peak p: candidate){
      result.add(new Peak(p.getMass()+Element.HYDROGEN_MASS,p.getIntensity(),0));
    }
    result.remove(result.size()-1);
    //Norm the distribution
    double intSum = 0.0;
    for (Peak p : result){
      intSum += p.getIntensity();	      
    }
    for (Peak p : result){
      p.setIntensity(p.getIntensity()/intSum);
    }
    return result;
  }
  
  private List<Peak> mix(List<Peak> candidate, List<Peak> shifted, double alpha){
    List<Peak> result = new ArrayList<Peak>(candidate.size());
    Iterator<Peak> sIter = shifted.iterator();
    for (Peak c : candidate){
      Peak s = sIter.next();
      double newMass = alpha*c.getMass()+(1-alpha)*s.getMass();
      double newInt = alpha*c.getIntensity()+(1-alpha)*s.getIntensity();
      result.add(new Peak(newMass, newInt ,0));
    }
    // These three lines assertion only
    double sum = 0.0;
    for (Peak p : result){ sum += p.getIntensity(); }
    if (sum < 1.0 - 1e8 || sum > 1.0 + 1e8) System.out.println("Error: sum="+sum);
    
    return result;
  }

}
