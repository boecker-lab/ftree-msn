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
  import java.io.FileReader;
  import java.io.IOException;
  import java.util.List;
  import java.util.ArrayList;
  import java.util.Collections;
  import java.util.regex.Pattern;


/**
 * @author Kerstin Scheubert
 */
public class MassErrorCheck {

    public static void main(String[] args) throws ElementNotFoundException, IOException {
      
      ElementTable et = ElementTable.generateCHNOPS();
      BufferedReader sFile = new BufferedReader(new FileReader(args[0]));
      readSpectra(sFile, et);
      sFile.close();
      
      

    } // end function
    
   
      static Compound readSpectra(BufferedReader file, ElementTable et) throws ElementNotFoundException, IOException{
      String name = "", str;
      Molecule formula = null;
      int mode = 0, collisionEnergy = 0;
      double tic = 0.0;
      List<Spectrum> spectra = new ArrayList<Spectrum>();
      List<Peak> peaks = null;
      Pattern peakPattern = Pattern.compile("\\d+\\.?\\d*\\s+\\d+\\.?\\d*e?\\d*.*");
      int i = 1;
      for (String line = file.readLine(); line!=null; line = file.readLine()){
        line = line.trim();
        if (line.contains("compound")) {
          name = line.substring(line.indexOf("compound")+9);
        }
        else if (line.contains("formula")) {
          str = line.substring(line.indexOf("formula")+8);
          try { formula = new Molecule(str, et); }
          catch (ElementNotFoundException e){
            throw new IOException("in line "+i+". "+e.getMessage());        
          }
          catch (Exception e) { 
            throw new IOException("in line "+i+". \""+str+"\" is not a valid sum formula");
          }
        }
        else if (line.contains("charge")){ 
          str = line.substring(line.indexOf("charge")+7);
          try { mode = Integer.parseInt(str); }
          catch (NumberFormatException e) { throw new IOException("in line "+i+". String \""+str+"\" could not be parsed to integer.");}
        }
        else if (line.contains("collision")){
          if (peaks != null && peaks.size() > 0){
            if (tic == 0.0) throw new IOException("in line "+i+". TIC of previous spectrum undefined, cannot calculate raw intensities.\nUse \"-nraw\" switch or correct input file.");
            spectra.add(new Spectrum(collisionEnergy, tic, peaks));
          }
          int begin = 11, end = line.length();
          if (line.contains("CE")) { begin = line.indexOf("CE")+3; }
          if (line.contains("eV")) { end = line.indexOf("eV")-1; }
          else if (line.contains("//")) { end = line.indexOf("//")-1; }
          str = line.substring(begin, end);
          str.trim();
          try { collisionEnergy = Integer.parseInt(str); }
          catch (NumberFormatException e) { throw new IOException("in line "+i+". String \""+str+"\" could not be parsed to integer.");}
          tic = 0.0;
          peaks = new ArrayList<Peak>();
        }
        else if (line.contains("tic")){
          str = line.substring(line.indexOf("tic")+4);
          try { tic = Double.parseDouble(str); }
          catch (NumberFormatException e) { throw new IOException("in line "+i+". String \""+str+"\" could not be parsed to double.");}       
        }
        else if (peakPattern.matcher(line).matches()){
          if (peaks == null) throw new IOException("in line "+i+". Peaks given before collision energy definition.");
          String[] values = line.split("\\s+");
          if (values.length < 2) throw new IOException("in line "+i+". Line seems to contain a peak, but does not.");
          try { peaks.add(new Peak(Double.parseDouble(values[0]), Double.parseDouble(values[1]), 0)); } // energy is 0 for now. Set later!
          catch (NumberFormatException e) {
            throw new IOException("in line "+i+". String \""+values[0]+"\" or \""+values[1]+"\" could not be parsed to double");
          }
          if (values.length > 2 && !values[2].equals("Parent") && !values[2].equals("Isotop") && !values[2].equals("Noise")) {
            Molecule desired=new Molecule(values[2], et);
            double isMass = peaks.get(peaks.size()-1).getMass()+Element.PROTON_MASS,
                   desiredMass=desired.getMass()+Element.HYDROGEN_MASS,
                   absError=Math.abs(isMass-desiredMass),
                   relError=absError/isMass*1e6;
            System.out.println(isMass+" "+desired+" "+desiredMass+" "+absError+" "+relError);
          }
        } else if (line.equals("")) {
          // do nothing
        } else {
          //System.out.println("Input file line "+i+": \""+line+"\" ignored.");
        }
        ++i;
      }
      
      if (tic == 0.0) throw new IOException("in line "+i+". TIC of previous spectrum undefined, cannot calculate raw intensities.\nUse \"-raw off\" switch or correct input file.");   
      spectra.add(new Spectrum(collisionEnergy, tic, peaks));
      
      Collections.sort(spectra);
      // Set the energy values
      i = 0;
      for (Spectrum s : spectra){
        for (Peak p : s.getPeaks()) {
          p.setEnergy(i);
        }
        ++i;
      }
      double mass = (formula != null)?formula.getMass():0.0;
      return new Compound(formula, name, mode, mass, 0.0, "", spectra); //focusedMass and instrument are unknown.
    }
    
}
