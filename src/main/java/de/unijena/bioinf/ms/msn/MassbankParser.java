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

import java.io.*;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

/**
 * @author Kerstin Scheubert
 */
public class MassbankParser{
	public static void main(String[] args) throws IOException {
	
		ElementTable et = new ElementTable();
		et.add(new Element("C", 12.0, 4));
		et.add(new Element("H", 1.007825, 1));
		et.add(new Element("N", 14.003074, 3));
		et.add(new Element("O", 15.994915, 2));
		et.add(new Element("P", 30.973762, 3));
		et.add(new Element("S", 31.972071, 2));
	
		
		BufferedReader reader;
		Molecule formula = null;
		String line, name = "", instrument = "", filename="6a_Methylprednisolone_", path="", outString="";; 
		String[] array;
		int mode = 0, collisionEnergy;
		double mass = 0.0, focusedMass = 0.0;
		
		TreeMap<String, Integer> map = new TreeMap<String, Integer>();
		boolean errorFlag = true; // starts with true, so no compound is added in the first loop.
		
		File dir = new File("/home/fr_86_bi/parser_fuer_dateien/mass_bank/");
		System.out.println(dir.listFiles());
		
		Map<String, String> outputfiles = new HashMap<String, String>();
		
		for(File f : dir.listFiles()) { 
			try {
			  reader = new BufferedReader(new FileReader(f));				
				line = reader.readLine();
		  	while (line != null && !line.contains("CH$NAME")){
		  	  line = reader.readLine();
		  	}
		  	name = line.substring(line.indexOf("CH$NAME")+9);
		  	
		  	if (!outputfiles.containsKey(name)){  
		  		outString="";
		  	  //If names of two spectra do not match, output the old compound and start a new one.
	  			/*if (!errorFlag) { // if there was no problem with the formula: add compound
	  				//compounds.add(new Compound(formula, name, mode, mass, focusedMass, instrument, spectra));
	  				oostream = new ObjectOutputStream(new FileOutputStream("../../compounds/"+numberToFixedLengthString(compoundsDone, 4)+".jos"));
	  				oostream.writeObject(new Compound(formula, name, mode, mass, focusedMass, instrument, spectra));
	  				oostream.close();
	  				map.put(name+" "+instrument+" "+mode, new Integer(compoundsDone));
	  				// The map stores the numbers of the compound output files.
		  			++compoundsDone;
		  		}*/
	  			//errorFlag = false;

			  	
			  	outString=outString+">compound "+name+"\n";
			  	
			  	while (line != null && !line.contains("CH$FORMULA")){
			  	  line = reader.readLine();
			  	}
					formula = new Molecule(line.substring(line.indexOf("CH$FORMULA")+12), et);
					outString += ">formula "+formula+"\n";
			  	while (line != null && !line.contains("CH$EXACT_MASS")){
			  	  line = reader.readLine();
			  	}
			  	mass = Double.parseDouble(line.substring(line.indexOf("CH$EXACT_MASS")+15));
			  	outString +=">parentmass "+ mass+ "Da \n";
			  	while (line != null && !line.contains("AC$INSTRUMENT")){
			  	  line = reader.readLine();
			  	}
					instrument = line.substring(line.indexOf("AC$INSTRUMENT")+15);
			  	while (line != null && !line.contains("AC$ANALYTICAL_CONDITION: MODE")){
			  	  line = reader.readLine();
			  	}				
					if (line.substring(line.indexOf("AC$ANALYTICAL_CONDITION: MODE")+30).contains("NEGATIVE")) mode = -1;
					else mode = 1;
					outString += ">charge "+ mode +"\n";
					outputfiles.put(name, outString);
				}
		  	
		  	while (line != null && !line.contains("AC$ANALYTICAL_CONDITION: COLLISION_ENERGY")){
		  	line = reader.readLine();
		  	}
		  	collisionEnergy = Integer.parseInt(line.substring(line.indexOf("AC$ANALYTICAL_CONDITION")+42, line.length()-3));
		  	outString = outputfiles.get(name);
		  	outString+="\n>collision "+collisionEnergy +"\n";
		  	  /*while (line != null && !line.contains("MS$FOCUSED_ION: PRECURSOR_M/Z")){
		  	  line = reader.readLine();
		  	}
		  	focusedMass = Double.parseDouble(line.substring(line.indexOf("MS$FOCUSED_ION")+30));*/
		  	while (line != null && !line.contains("PK$PEAK")){
		  	line = reader.readLine();
		  	}
				line = reader.readLine();
				//peaks = new Vector<Peak>();
				while (line != null && !line.contains("//")){
					array = line.split(" ");
					// array[2] is mass, array[3] abs. intensity, array[4] rel. intensity.
					// spectra.size shows how many spectra had a lower energy than the spectrum this peak belongs to.
					//peaks.add(new Peak(Double.parseDouble(array[2]), Double.parseDouble(array[3]), Double.parseDouble(array[4]), spectra.size()));
					line = reader.readLine();
					outString+= array[2] + " " +array[3] + "\n";
				}
				
				outputfiles.put(name, outString);
				
				//spectra.add(new Spectrum(collisionEnergy, peaks));
				
			}
			catch (ElementNotFoundException e) {
				//System.out.println("Spectrum No."+i+": "+e.getMessage());
				errorFlag = true;
			}
			
		}
		
		System.out.println("Reading done");
		/*oostream = new ObjectOutputStream(new FileOutputStream("/home/fr_86_bi/parser_fuer_dateien/massbank/hill/index.jos"));
		oostream.writeObject(map);
		oostream.close();
		System.out.println("Map written");*/		
		for (String outName : outputfiles.keySet()){
		  path =  "/home/fr_86_bi/parser_fuer_dateien/massbank/hill/"+outName+".ms";
		  try {
	      BufferedWriter out = new BufferedWriter(new FileWriter(path));
	      out.write(outputfiles.get(outName)); System.err.println("result written to file: "+path);
	      out.close();
	    } catch (IOException e) {
	    }
		}
		
				
	}
	
}
