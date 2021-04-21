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
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * @author Kerstin Scheubert
 */
public class MsToMbConverter {

  public static void main(String[] args) throws Exception{
    File dir = new File(args[0]);
    Parameters param = new Parameters();
    param.raw = false;

    BufferedReader tableFile = new BufferedReader(new FileReader(args[1]));
    // Read table
    tableFile.readLine(); //discard header
    Map<String, MM48Compound> table = new HashMap<String, MM48Compound>();

    //this is for mm48
    //header should be: No. Name  Formel  Ladung  Compound ID Retentionszeit [s]  M [M+H]+  [M]+  [M+Na]+ Fragment
    for (String line = tableFile.readLine(); line!=null; line = tableFile.readLine()){
      String[] values = line.split("\\t");
      table.put(values[1], new MM48Compound(values[1], values[2], Double.parseDouble(values[5]), Double.parseDouble(values[6]), values[4]));
    }

    // header should be Name;CID;Formula;Weight;
    /*for (String line = tableFile.readLine(); line!=null; line = tableFile.readLine()){
      String[] values = line.split(";");
      table.put(values[0], new MM48Compound(values[0], values[2], 0.0, 0.0, values[1]));
    } */   

    for (String f: dir.list()) {
      if (f.endsWith(".ms")) {
        BufferedReader file = new BufferedReader(new FileReader(args[0]+"/"+f));
        String o = f.substring(0, f.indexOf('.')) + ".mb";
        BufferedWriter out = new BufferedWriter(new FileWriter(args[0]+"/"+o));
        Compound c = MS2Analyzer.readSpectra(file, param);
        String key = c.getName().trim().replaceAll("_", " ");

        file.close();

        List<Peak> peaks = c.mergePeaks(param);
        out.write("# Sample: ");
        out.write(c.getName());
        out.newLine();
        out.write("# eV: 0\n");
        out.write("# Pseudospectrum: 0\n");
        out.write("# Retentiontime: 0\n");
        out.write("# Parent Mass: ");
        out.write(Double.toString(c.getMass()));
        out.newLine();
        out.write("# PubChem ID: ");
        try{
          out.write(table.get(key).cid);
        } catch (RuntimeException e) {                    
          System.err.println(f.toString());
          e.printStackTrace();
        }
        out.newLine();
        for (Peak p : peaks) {
          out.write(p.getMass()+" "+p.getIntensity());
          out.newLine();
        }
        out.close();
      }      

    }    

  }

}
