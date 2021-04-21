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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;


/**
 * @author Kerstin Scheubert
 */
public class SvatosOptimizer {

  /**
   * @param args
   */
  public static void main(String[] args) throws IOException  {
    Parameters param = new Parameters();
    param.intenseMerge = true;
    File dir = new File(args[0]);
    
    for (String f: dir.list()) {
      if (f.endsWith(".ms")) {
        BufferedReader file = new BufferedReader(new FileReader(args[0]+"/"+f));
        Compound c = MS2Analyzer.readSpectra(file, param);
        double minCE = Double.POSITIVE_INFINITY;
        for (Spectrum s : c.getSpectra()){
          if (s.getCollisionEnergy() < minCE){
            minCE = s.getCollisionEnergy();
          }
        }
        double maxTIC = 0;
        int argmaxTIC = 0, i=0;
        for (Spectrum s: c.getSpectra()){
          if (s.getCollisionEnergy() == minCE){
            if (s.getTic() > maxTIC){
              maxTIC = s.getTic();
              argmaxTIC = i;
            }
          }
          ++i;
        }
        
        List<Spectrum> bestSpectra = new ArrayList<Spectrum>();
        Compound result = new Compound(c.getFormula(), c.getName(), c.getCharge(), c.getMass(), c.getFocusedMass(), c.getInstrument(), bestSpectra);
        boolean CeIncrease = true;
        for (int j = argmaxTIC; CeIncrease && j < c.getSpectra().size(); ++j){
          bestSpectra.add(c.getSpectra().get(j));
          if (j < c.getSpectra().size()-1 && c.getSpectra().get(j).getCollisionEnergy() > c.getSpectra().get(j+1).getCollisionEnergy()){
            CeIncrease = false;
          }
        }
        
        for (Spectrum s : bestSpectra){
          s.mergePeaks(param);
        }
        BufferedWriter out = new BufferedWriter(new FileWriter(args[1]+"/"+f));
        out.write(">compound "+result.getName());
        out.write("\n>formula "+result.getFormula().toString());
        out.write("\n>parentmass "+result.getFocusedMass());
        out.write("\n>charge 1\n\n");
        for (Spectrum s : result.getSpectra()){
          out.write("#>retention "+s.getRT());
          out.write("\n>collision "+s.getCollisionEnergy());
          out.write("\n>tic "+s.getTic());
          out.newLine();
          for (Peak p : s.getPeaks()){
            out.write(p.getMass()+" "+p.getIntensity());
            out.newLine();
          }
          out.newLine();          
        }
        out.close();
        file.close();
      } // end if ms-file
    } // end for all files in dir

  }

}
