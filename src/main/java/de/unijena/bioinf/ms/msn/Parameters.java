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

import java.io.Writer;
import java.io.IOException;

/**
 * @author Kerstin Scheubert
 */
public class Parameters{
	public Compound input;
	public ElementTable et;
	public double errorPpm;
	public double absError;
	public double parentPeakError;
	public double intError;
	public double intOffset;
  public double errorIncrease; // multiplicative increase of the error that should be applied at zero intensity  
	public double precision;
	public boolean raw;
	public boolean relative;
	public double mergeThreshold;
	public boolean intenseMerge;
	public double intensityCutoff;
	public boolean scan;
	public boolean isotopes;
  public double isobonus;
  public boolean removeiso;
  public boolean parentiso;
	public boolean gheuristic;
	public boolean tdheuristic;
  public int peaksToScore;
  public double[] closureScore;
	public double massDeviationPenalty;
	public String neutralLossList;
	public int nlCombinations;
	public boolean hcRatioScoring;
	public boolean heteroScoring;
	public boolean DBEDependentScoring;
	public double hcAverage;
	public double hcDev;
	public double heteroAverage;
	public double heteroDev;
	public double upToXDa; // perform DBE-Scoring only up to this mass
	public double DBEAverage;
	public double DBEStdDev;
	public String graphfile;
	public boolean writeGraph;
	public boolean writeGraphWithTransitiveClosure;
	public String outfile;
	public int detailedOutput;
	public boolean success;
	  public String radicalLossList;
	  public String forbiddenLossList;
	public boolean attachRemaining;	
	public boolean mostIntensivePeaks;
	public boolean noTransitiveClosureScore;
	public boolean runningTimeAnalysis;
	public String runningTimeAnalysisFile;
	public boolean showTransitiveClosure;
	public boolean integerprog;
	
	// Fill the values with defaults, the required ones with null
	public Parameters(){
		input=null;
		et = ElementTable.generateCHNOPS();
		absError = 1.0;
		parentPeakError=0.5;
		errorPpm = 20.0;
		intError = 3.0;
		intOffset = 0.02;
    errorIncrease = 1.0;
		precision = 1e-5;
		raw = true;
		relative = false;
		mergeThreshold = 0.1;
		intenseMerge = false;
		attachRemaining = false;
		mostIntensivePeaks=false;
		noTransitiveClosureScore=false;		
		intensityCutoff = 0;
		scan = false;
		isotopes = false;
    isobonus = 0;
    parentiso=false;
    removeiso=false;
		gheuristic = false;
		integerprog= false;
		tdheuristic = false;
    peaksToScore = 10;
    closureScore=new double[3];
    closureScore[0]=3;
    closureScore[1]=-0.5;
    closureScore[2]=-0.5;
		massDeviationPenalty = 3.0;
		neutralLossList = "C2H2 C2H4 C2H4O2 C2H5O4P C3H2O3 C3H4O4 C3H6O2 C3H9N C4H8 C5H8 C5H8O4 C6H10O4 C6H10O5 C6H6 C6H8O6 CH2O CH2O2 CH3 CH3N CH4 CH4N2O CH4O CH5N CHNO CO CO2 H2 H2O H2S H2SO4 HPO3 N2 NH3 S SO2 SO3";
	    radicalLossList = "H OH CH3 CH3O C3H7 C4H9 C6H5O O";
	    forbiddenLossList = "H2 C2O C4O C3H2 C5H2 C7H2";
	    //neutralLossList = "";
    nlCombinations = 3;
		hcRatioScoring = false;
		DBEDependentScoring = false;
		upToXDa = 200;
		DBEAverage = 6.151312; //new from knapsack: 9 old from KEGG: 8.143
		DBEStdDev = 4.541604;
		heteroScoring = false;
		hcAverage = 1.435877; // Gut: 1.3 0.8 Aracyc: 1.6 0.5 Knapsack: 1.3 0.4 Alt: 1.2 0.8
		hcDev = 0.4960778;
		heteroAverage = 0.5886335; // new: 0 1 old: 0.2 0.5
		heteroDev = 0.5550574;		
		detailedOutput = 10;
		graphfile = "MS2Graph";
		writeGraph = false;
		writeGraphWithTransitiveClosure=false;
		outfile = "MS2Analyzer.out";
		runningTimeAnalysis=false;
		runningTimeAnalysisFile="runningtime.txt";
		success = false;
		showTransitiveClosure=false;
	}

  public void write(Writer out) throws IOException{
    StringBuffer b = new StringBuffer();
    b.append("Alphabet: ");
    b.append(et);
    b.append("\nErrors: ");
    b.append(absError);
    b.append(" mDa, ");
    b.append(errorPpm);
    b.append(" ppm, Intensity: ");
    b.append(intError);
    b.append(" %, Increase: ");
    b.append(errorIncrease);
    b.append("\n Intensity conversion:");
    if (raw){
      b.append(" pseudo raw");
    } else if (relative){
      b.append(" relative/ranked");
    } else {
      b.append(" none");
    }
    b.append(" Merge threshold: ");
    b.append(mergeThreshold);
    if (intenseMerge){
      b.append(" (retained most intense peak)");
    }
    b.append("\n");
    if (isotopes){
      b.append("Used Isotopes ");
    }
    if(gheuristic){
      b.append("Used greedy heuristic ");
    }
    if(tdheuristic){
      b.append("Used top-down heuristic ");
    }
    b.append("\nPeaks considered in exact calculation: ");
    b.append(peaksToScore);
    b.append("\nScoring: Mass dev. ");
    b.append(massDeviationPenalty);
    b.append(" Neutral loss comb. ");
    b.append(nlCombinations);
    b.append("\nNeutral loss list: ");
    b.append(neutralLossList);
    if (hcRatioScoring){
      b.append("\nHC Average: ");
      b.append(hcAverage);
      b.append(" Dev: ");
      b.append(hcDev);      
    }
    if (heteroScoring){
      b.append("\nHeteroRatio Average: ");
      b.append(heteroAverage);
      b.append(" Dev: ");
      b.append(heteroDev);      
    }
    if (DBEDependentScoring){
      b.append("\nDBE Average: ");
      b.append(DBEAverage);
      b.append(" Dev: ");
      b.append(DBEStdDev);
      b.append("If fragment was not heavier than: ");
      b.append(upToXDa);
    }
    out.write(b.toString());
  }
  
  public double getErrorForMass(double mass){
  	return Math.max(errorPpm*1e-6*mass, absError*1e-3);
  }
}
