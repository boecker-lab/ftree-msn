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
import java.util.ArrayList;

/**
 * @author Kerstin Scheubert
 */
public class IsotopeDistribution {
    
    public static final int ISOTOPE_DISTRIBUTION_SIZE = 10;
    private Integer nominalMass;
    
    private List<Peak> peaks = new ArrayList<Peak>(ISOTOPE_DISTRIBUTION_SIZE);
	
    public IsotopeDistribution(Integer nominalMass, List<Peak> peaks) {
        this.nominalMass = nominalMass;
        this.peaks = peaks;
    }
    
    public IsotopeDistribution(Integer nominalMass, double[][] peaks) {
		this.nominalMass = nominalMass;
        array2list(peaks);
	  }

	public IsotopeDistribution() {
		this(new Integer(0), new ArrayList<Peak>(4));
	}

	public void addIsotope(Peak p){
	  peaks.add(p);
	}
	
	private void array2list(double[][] peaksArray) {
		if (peaksArray != null) {
			for (int i = 0; i < peaksArray[0].length && i < IsotopeDistribution.ISOTOPE_DISTRIBUTION_SIZE; ++i) {
				peaks.add(new Peak(peaksArray[0][i], peaksArray[1][i], 0));
			}			
		}						
	}

	private double[][] foldquick(final double[][] peaks, int power) {
    	if (power == 1) {
    		return peaks;
    	}
    	if (peaks == null) {
    		return null;
    	}
    	
    	int binarySize = (int)Math.ceil(Math.log((double)power) / Math.log(2)) + 1;
    	int[] binary = new int[binarySize];
    	int i = 0;
    	while (power > 0) {
    		binary[i++] = power & 1;
    		power >>= 1;
    	}
    	if (binary.length == 0) {
    		return peaks;
    	}
    	    	
    	double[][] powerPeaks = new double[2][IsotopeDistribution.ISOTOPE_DISTRIBUTION_SIZE];
    	System.arraycopy(peaks[0], 0, powerPeaks[0], 0, IsotopeDistribution.ISOTOPE_DISTRIBUTION_SIZE);
    	System.arraycopy(peaks[1], 0, powerPeaks[1], 0, IsotopeDistribution.ISOTOPE_DISTRIBUTION_SIZE);
    	
    	double[][] resultPeaks = null;
    	if (binary[0] == 1) {
    		resultPeaks = new double[2][IsotopeDistribution.ISOTOPE_DISTRIBUTION_SIZE];
        	System.arraycopy(peaks[0], 0, resultPeaks[0], 0, IsotopeDistribution.ISOTOPE_DISTRIBUTION_SIZE);
        	System.arraycopy(peaks[1], 0, resultPeaks[1], 0, IsotopeDistribution.ISOTOPE_DISTRIBUTION_SIZE);
    	}
    	for (int j = 1; j < binary.length; j++) {
			powerPeaks = foldquick(powerPeaks, powerPeaks);
			if (binary[j] == 1) {
				resultPeaks = foldquick(resultPeaks, powerPeaks);
			}
		}
    	return resultPeaks;
    }
    
    private double[][] foldquick(final double[][] peaks1, final double[][] peaks2) {
    	
    	if (peaks1 == null) {
    		return peaks2;
    	}
    	
    	if (peaks2 == null) {
    		return peaks1;
    	}
    	
    	double[][] result_peaks = new double[2][IsotopeDistribution.ISOTOPE_DISTRIBUTION_SIZE];
    	
        double abundances_sum, masses_mult_abundances_sum;
        
        for (int i = 0; i < IsotopeDistribution.ISOTOPE_DISTRIBUTION_SIZE; ++i) {
            abundances_sum = 0.0;
            masses_mult_abundances_sum = 0.0;
            int i1 = 0;
            int i2 = i + 1;
            
            while (--i2 >= 0) {
            	abundances_sum += peaks1[1][i1] * peaks2[1][i2];
            	masses_mult_abundances_sum += peaks1[1][i1] * peaks2[1][i2] * (peaks1[0][i1] + peaks2[0][i2]);
            	++i1;
            }
            
            double mass = (abundances_sum != 0) ? 
                    masses_mult_abundances_sum / abundances_sum: 0.0;
            result_peaks[0][i] = mass;
            result_peaks[1][i] = abundances_sum;
        }
    	
    	return result_peaks;
    }
    
    public IsotopeDistribution fold(int power) {
        if (power <= 1) {
            return this;
        }
        if (isEmpty()) {
        	return this;
        }
        
        double[][] fastPeaks = new double[2][IsotopeDistribution.ISOTOPE_DISTRIBUTION_SIZE];

        // copies first everything from peaks
        int i = 0;
        for (; i < IsotopeDistribution.ISOTOPE_DISTRIBUTION_SIZE && i < peaks.size(); ++i) {
			fastPeaks[0][i] = peaks.get(i).getMass();
			fastPeaks[1][i] = peaks.get(i).getIntensity();
		}
        // fills the rest (if necessary) until default size
        for (; i < IsotopeDistribution.ISOTOPE_DISTRIBUTION_SIZE; ++i) {
			fastPeaks[0][i] = 0.0;
			fastPeaks[1][i] = 0.0;
		}
        
        double[][] resultPeaks = foldquick(fastPeaks, power);
                
        return new IsotopeDistribution(nominalMass*power, resultPeaks);
    }
	
	public IsotopeDistribution fold(IsotopeDistribution distribution) {
        
        if (distribution.isEmpty()) {
            return this;
        }
        if (isEmpty()) {
            return distribution;
        }
        
		double[][] fastPeaks1 = new double[2][IsotopeDistribution.ISOTOPE_DISTRIBUTION_SIZE];
        
        // copies first everything from peaks
		int i = 0;
        for (; i < IsotopeDistribution.ISOTOPE_DISTRIBUTION_SIZE && i < peaks.size(); ++i) {
			fastPeaks1[0][i] = peaks.get(i).getMass();
			fastPeaks1[1][i] = peaks.get(i).getIntensity();
		}
        // fills the rest (if necessary) until default size
        for (; i < IsotopeDistribution.ISOTOPE_DISTRIBUTION_SIZE; ++i) {
			fastPeaks1[0][i] = 0.0;
			fastPeaks1[1][i] = 0.0;
		}

		double[][] fastPeaks2 = new double[2][IsotopeDistribution.ISOTOPE_DISTRIBUTION_SIZE];

        // copies first everything from peaks
		i = 0;
        for (; i < IsotopeDistribution.ISOTOPE_DISTRIBUTION_SIZE && i < distribution.size(); ++i) {
			fastPeaks2[0][i] = distribution.getPeak(i).getMass();
			fastPeaks2[1][i] = distribution.getPeak(i).getIntensity();
		}

        // fills the rest (if necessary) until default size        
        for (; i < IsotopeDistribution.ISOTOPE_DISTRIBUTION_SIZE; ++i) {
			fastPeaks2[0][i] = 0.0;
			fastPeaks2[1][i] = 0.0;
		}
        
        
		double[][] resultPeaks = foldquick(fastPeaks1, fastPeaks2);
		
		return new IsotopeDistribution(nominalMass + distribution.getNominalMass(), resultPeaks);
	}
	
	
    /**
     * @return Returns the peaks.
     */
    public List<Peak> getPeaks() {
        return peaks;
    }
    
    /**
     * @param peaks The peaks to set.
     */
    public void setPeaks(List<Peak> peaks) {
        this.peaks = peaks;
    }
        
    /**
     * @return Returns the nominalMass.
     */
    public Integer getNominalMass() {
        return nominalMass;
    }
    /**
     * @param nominalMass The nominalMass to set.
     */
    public void setNominalMass(Integer nominalMass) {
        this.nominalMass = nominalMass;
    }
    
    public boolean isEmpty() {
        return peaks.isEmpty();
    }

	public int size() {
		return peaks.size();
	}

	public Peak getPeak(int i) {
		return peaks.get(i);
	}
	
}
