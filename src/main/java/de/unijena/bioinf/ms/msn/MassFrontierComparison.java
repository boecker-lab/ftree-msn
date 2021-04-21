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
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;


/**
 * @author Kerstin Scheubert
 */
public class MassFrontierComparison {

  public static void main(String[] args) {
    
    ElementTable et = new ElementTable();
    et.add(new Element("C", 12.0, 4));
    et.add(new Element("H", 1.007825, 1));
    et.add(new Element("N", 14.003074, 3));
    et.add(new Element("O", 15.994915, 2));
    et.add(new Element("P", 30.973762, 3));
    et.add(new Element("S", 31.972071, 2));
    
    boolean createfile=false;
    boolean readfile=true;
    String path = "/home/m3rafl/projects/metabo_tandem_ms/data/hill/MassFrontierResults/allCE.csv";
    String text=""; 
    String line, name="";
      double mass;
      BufferedReader reader;
      
      Map<String,List<Double>> map = new HashMap<String, List<Double>>();
      if(createfile){
        for (int i=1; i<=5; i++){
          try
          {
            reader = new BufferedReader(new FileReader("/home/m3rafl/projects/metabo_tandem_ms/data/hill/MassFrontierResults/"+i+"0eV.csv"));   
            System.out.println("reading file: /home/m3rafl/projects/metabo_tandem_ms/data/hill/MassFrontierResults/"+i+"0eV.csv");
            line = reader.readLine();
            while(line.contains(",")){
              List<Double> masses = new ArrayList<Double>();
              int endIndex=line.indexOf(",");
                name = line.substring(1,endIndex-1);
                System.out.println(name);
                if(map.containsKey(name)){
                  masses=map.get(name);
                  System.out.println("map contains the molecule");
                }
                    
                int beginIndex;
                while(!(line.substring(endIndex+1).startsWith(","))&& endIndex!=-1){
                  beginIndex=endIndex;
                  endIndex=line.indexOf(",",beginIndex+1);
                  if (endIndex==-1){
                    mass=Double.parseDouble(line.substring(beginIndex+1));
                  }else{
                    mass=Double.parseDouble(line.substring(beginIndex+1,endIndex));
                  }
                  
                  System.out.println(mass);
                  if(!(masses.contains(mass))){
                    masses.add(mass);
                  }
                  
                }
                
                System.out.println(masses);
                map.put(name, masses);
                line = reader.readLine();

            }

          }
          catch(Exception eve)
          {
          }
        }
        
        
        
        for (Iterator iter = map.keySet().iterator(); iter.hasNext();) {
          Object tmpKey = iter.next();
          Object tmpValue = map.get(tmpKey);

          text+=tmpKey.toString();
          text+=tmpValue.toString();
          text+="\n";
          System.out.println("Key: " + tmpKey.toString());
          //System.out.println("Value: " + tmpValue.toString());
          //System.out.println();
          }
        
        //write to file
        try {
              BufferedWriter out = new BufferedWriter(new FileWriter(path));
              out.write(text); System.err.println("result written to file: "+path);
              out.close();
          } catch (IOException e) {
          }
      }
      
       
       if (readfile){
         try{
           reader = new BufferedReader(new FileReader("/home/m3rafl/projects/metabo_tandem_ms/data/hill/MassFrontierResults/allCE.csv"));   
           line = reader.readLine();
           while(line.contains("[")){
             List<Double> masses = new ArrayList<Double>();
             int endIndex=line.indexOf("[");
               name = line.substring(0,endIndex);
               System.out.println(name);
                   
               int beginIndex;
               while(!(line.substring(endIndex).startsWith("]"))&& endIndex!=-1 && !(line.substring(endIndex+1).startsWith("]"))){
                 beginIndex=endIndex;
                 endIndex=line.indexOf(",",beginIndex+1);
                 if (endIndex==-1){
                   endIndex=line.indexOf("]",beginIndex+1);
                   mass=Double.parseDouble(line.substring(beginIndex+1, endIndex));
                 }else{
                   mass=Double.parseDouble(line.substring(beginIndex+1,endIndex));
                 }
                 
                 System.out.println(mass);
                 if(!(masses.contains(mass))){
                   masses.add(mass);
                 }
                 
               }
               
               System.out.println(masses);
               map.put(name, masses);
               line = reader.readLine();
           }
         }
         catch(Exception eve)
         {
         }
       }
        
      
        File dir = new File("/home/m3rafl/projects/metabo_tandem_ms/data/hill/treesforMF");
      //System.out.println(dir.listFiles());
      
      
        String pathExact = "/home/m3rafl/projects/metabo_tandem_ms/data/hill/exact.txt";
        String pathUnexact = "/home/m3rafl/projects/metabo_tandem_ms/data/hill/unexact.txt";
        String pathNotfound = "/home/m3rafl/projects/metabo_tandem_ms/data/hill/notfound.txt";
        
        String textExact="name, MF results, tree, nr of results, results \n";
        String textUnexact="name, MF results, tree, nr of results, results \n";
        String textNotfound="name, MF results, tree, nr of results, results \n";
        
        for(File f : dir.listFiles()) { 
        try {
          reader = new BufferedReader(new FileReader(f));
          line = reader.readLine();
          double moleculeMass;
          Molecule formula = null;
          String molecule="";
          List<Double> massList = new ArrayList<Double>();
          //Ueberspringen des Anfangs bis zu den Molekuelen
          while (line != null && !line.contains("label")){
              line = reader.readLine();
          }
          //Name einlesen
          name=line.substring(line.indexOf("for")+4, line.indexOf("correct")-1);
          if(name.contains("_")){
            name.indexOf("_");
            String[] newName = name.split("_");
            name=newName[0];
            for (int i=1; i<newName.length; i++){
              name+=" "+newName[i];
            }
          }
          System.out.println(name);
          while (line != null && !line.contains("[")){
              line = reader.readLine();
          }
          
          //Liste der Massen erstellen
          while(line!=null && line.contains("label")){
            if(!(line.contains("->"))){
              molecule=line.substring(0,line.indexOf("["));
              formula=new Molecule(molecule,et);
              moleculeMass= formula.getMass();
              //System.out.println(moleculeMass);
              if(!(massList.contains(moleculeMass))){
                massList.add(moleculeMass);
              }
              
            }
            line=reader.readLine();
          }
          
          System.out.println("massList: "+massList);
          
          List<Double> compareList = map.get(name);
          System.out.println("compareList: "+compareList);
          
          List<Double> exactList = new ArrayList<Double>();
          List<Double> unexactList = new ArrayList<Double>();
          List<Double> notfoundList = new ArrayList<Double>();
          
          
          
          textExact+=name+", "+compareList.size()+", "+massList.size()+", ";
          textUnexact+=name+", "+compareList.size()+", "+massList.size()+", ";
          textNotfound+=name+", "+compareList.size()+", "+massList.size()+", ";
          
          for (Iterator iterCompareList =compareList.iterator(); iterCompareList.hasNext(); ){
              double compareMass = (Double) iterCompareList.next();
              //System.out.println(compareMass);
              boolean isExact = false, isUnexact = false;
              for (Iterator iterMassList = massList.iterator(); iterMassList.hasNext();){
                  
                  double exactMass=(Double) iterMassList.next();
                  exactMass = exactMass + Element.PROTON_MASS;
                  double exactMass2 = Math.round((exactMass+Element.HYDROGEN_MASS) * 10000.)/10000.,
                  exactMass3 = Math.round((exactMass+2*Element.HYDROGEN_MASS) * 10000.)/10000.;
                  exactMass = Math.round(exactMass * 10000.)/10000.;
                  //System.out.println(exactMass);
                  if(compareMass==exactMass || compareMass==exactMass2 || compareMass == exactMass3){
                    //exakte Treffer
                    isExact = true;                   
                  }
                  else if(Math.abs(compareMass-exactMass) < 0.5 || Math.abs(compareMass-exactMass2) < 0.5 || Math.abs(compareMass-exactMass3) < 0.5){
                    //unexakte Treffer
                    isUnexact = true;
                  }
              }
              
              if (isExact){
                exactList.add(compareMass);
              } else if (isUnexact){
                unexactList.add(compareMass);
              } else {                
                notfoundList.add(compareMass);
              }
          }
          
          textExact+=exactList.size()+", " +unexactList.size()+", "+notfoundList.size()+"\n";
          textUnexact+=unexactList.size()+", "+unexactList+"\n";
          textNotfound+=notfoundList.size()+", "+notfoundList+"\n";
          
        }
        catch(Exception eve){
        }
        
        
      }
        try {
          BufferedWriter out = new BufferedWriter(new FileWriter(pathExact));
          out.write(textExact); 
          //System.err.println("result written to file: "+pathExact);
          out.close();
        } catch (IOException e) {
        }
        try {
          BufferedWriter out = new BufferedWriter(new FileWriter(pathUnexact));
          out.write(textUnexact); 
          //System.err.println("result written to file: "+pathUnexact);
          out.close();
        } catch (IOException e) {
        }
        try {
          BufferedWriter out = new BufferedWriter(new FileWriter(pathNotfound));
          out.write(textNotfound); 
          //System.err.println("result written to file: "+pathNotfound);
          out.close();
        } catch (IOException e) {
        }
    }

}
