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
import java.util.regex.Pattern;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Kerstin Scheubert
 */
public class AnnotateMsFiles {

	/**
	 * @param args 1. File with MM48 contents 2. Input dir 3. Output dir
	 */
	public static void main(String[] args) throws IOException {
		BufferedReader tableFile = new BufferedReader(new FileReader(args[0]));
		File dir = new File(args[1]);
		// Read table
		tableFile.readLine(); //discard header
		//header should be: No.	Name	Formel	Ladung	Compound ID	Retentionszeit [s]	M	[M+H]+	[M]+	[M+Na]+	Fragment
		List<MM48Compound> table = new ArrayList<MM48Compound>();
		for (String line = tableFile.readLine(); line!=null; line = tableFile.readLine()){
			String[] values = line.split("\\t");
			table.add(new MM48Compound(values[1], values[2], Double.parseDouble(values[5]), Double.parseDouble(values[6]), values[4]));
		}
		
		for (String fileStr : dir.list()){
			if (fileStr.endsWith("ms")){
				File f = new File(args[1]+"/"+fileStr);
				BufferedReader file = new BufferedReader(new FileReader(f));
				StringBuffer rest = new StringBuffer();
        double focusedMass = 0.0, retTime = 0.0;
        boolean parentPeakSeen = false, peaksSeen = false;
    		Pattern peakPattern = Pattern.compile("\\d+\\.?\\d*\\s+\\d+\\.?\\d*e?\\d*");        
				for (String line = file.readLine(); line!=null; line = file.readLine()){
					line = line.trim();
					if (line.contains("compound")) {
						// Do nothing
					}
		      else if (line.contains("parentmass")){ 
		        String str = null;
		      	if (line.contains("Da")){
		          str = line.substring(line.indexOf("parentmass")+11, line.indexOf("Da")-1);
		        } else {
		          str = line.substring(line.indexOf("parentmass")+11);          
		        }
		        try { focusedMass = Double.parseDouble(str); }
		        catch (NumberFormatException e) { throw new IOException("in file "+fileStr+". String \""+str+"\" could not be parsed to double.");}
		      	rest.append(line);
		      	rest.append("\n");
		      } else if (line.contains("retention")){ 
			        String str = null;
			      	if (line.contains("s")){
			          str = line.substring(line.indexOf("parentmass")+12, line.indexOf("s")-1);
			        } else {
			          str = line.substring(line.indexOf("parentmass")+12);          
			        }
			        try { retTime = Double.parseDouble(str); }
			        catch (NumberFormatException e) { throw new IOException("in file "+fileStr+". String \""+str+"\" could not be parsed to double.");}
			      	rest.append(line);
			      	rest.append("\n");

		      } else {
		      	rest.append(line);
		      	rest.append("\n");
		      	if (peakPattern.matcher(line).matches()){
		  				peaksSeen = true;
		      		String[] values = line.split("\\s+");
		  				if (values.length < 2) throw new IOException("in file "+fileStr+". Line seems to contain a peak, but does not.");
		  				double mass = Double.parseDouble(values[0]);
		      		if (Math.abs(mass-focusedMass)< 1.0){
		      			parentPeakSeen = true;
		      		}
		      	}
		      }
				}
				file.close();
				
				if (peaksSeen) {
					String ppString = parentPeakSeen ? "" : "npp";
					MM48Compound best = null;
				  final double protonMass = 1.007276;
					// determine name and formula from table
					for (MM48Compound candidate : table){
						double chargeMass = candidate.formula.contains("+")?0.0:protonMass; 
						if (Math.abs(candidate.mass-focusedMass+chargeMass)<1){
							if(best == null || best.calcDiff(retTime, focusedMass) < candidate.calcDiff(retTime, focusedMass)){
								best = candidate;
							}
						}
					}
					
					
					if (best != null) {
						File outfile = new File(args[2] + "/" + best.name + ppString + ".ms");
						int i = 2;
						while (outfile.exists()) {
							outfile = new File(args[2] + "/" + best.name + i + ppString + ".ms");
							++i;
						}
						BufferedWriter out = new BufferedWriter(new FileWriter(outfile));
						out.write(">compound ");
						out.write(best.name);
						out.newLine();
						out.write(">formula ");
						out.write(best.formula);
						out.newLine();
						out.write(rest.toString());
						out.close();
					} else {
						System.out.println("No hit for "+fileStr);					
					}
				}
				
			}
		}
	}

}

class MM48Compound{
	public double retTime;
	public double mass;
	public String name;
	public String formula;
  public String cid;
	
	public MM48Compound(String name, String formula, double retTime, double mass, String cid) {
		super();
		this.name = name;
		this.formula = formula;
		this.retTime = retTime;
		this.mass = mass;
    this.cid = cid;
	}
	
	public double calcDiff(double retTime, double mass){
		return Math.abs(this.retTime-retTime);//+5*Math.abs(this.mass-mass);
	}
	
	public String toString(){
		return name+" "+formula+" "+retTime+" "+mass+"\n";
	}
}