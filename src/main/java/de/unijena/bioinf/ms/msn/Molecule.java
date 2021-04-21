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

import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Arrays;

/**
 * @author Kerstin Scheubert
 */
public class Molecule implements java.io.Serializable, Comparable<Molecule> {

  private static final long serialVersionUID = -2142593192792902819L;
  private Abundancy[] contents;
	private int position;
	private double mass;
	private ElementTable elementTable;
	private double dbe;
	private double score;
  private double isotopeScore;
	private Vertex vertex;
	
	public Molecule(Molecule another){
		elementTable = another.elementTable;
		position = another.position;
		mass = another.mass;
		dbe = another.dbe;
		score = another.score;
		isotopeScore = another.isotopeScore;
		vertex = null;
		contents = new Abundancy[elementTable.size()];
		for (int i = 0; i < position; ++i){
			contents[i] = new Abundancy(another.contents[i].name, another.contents[i].amount);
		}
	}
	
	public Molecule(ElementTable et){
		mass = 0.0;
		elementTable = et;
		contents = new Abundancy[et.size()];
		position = 0;
		score = 0.0;
    isotopeScore = 0;
		dbe = 1.0;
		vertex = null;
	}
	
	public Molecule(ElementTable et, int[] compomer, Weights elOrder) throws ElementNotFoundException {
		this(et);
		if (compomer.length != elOrder.length()) System.out.println("Length of Compomer and Weights have to be equal!");
		for (int i = 0; i < compomer.length; ++i){
			if (compomer[i] > 0){
				mass += et.get(elOrder.getName(i)).getMass()*compomer[i];
				contents[position] = new Abundancy(elOrder.getName(i), compomer[i]);
				++position;
			}
		}
		Arrays.sort(contents, 0, position);
		dbe = calcDBE();
	}
	
	public Molecule(String formula, ElementTable et) throws ElementNotFoundException {
		this(et);
		StringBuffer element = new StringBuffer(), number = new StringBuffer();
		Element el;
		int amount;
		for (int i=0; i < formula.length(); ++i){
			char c = formula.charAt(i);
			if (Character.isUpperCase(c)){
				// if an element has been read
				if (element.length() != 0) {
					// lookup element
					el = et.get(element.toString());
					// if a number has been read: parse amount, enter in map and update mass.
					if (number.length() != 0) {
						amount = Integer.parseInt(number.toString());
						contents[position] = new Abundancy(element.toString(), amount);
						++position;
						mass += amount*el.getMass();
					}	else { // Insert one element
						contents[position] = new Abundancy(element.toString(), 1);
						++position;
						mass += el.getMass();
					}
				}
				element.delete(0, element.length());
				number.delete(0, number.length());
				element.append(c);
			} else if (Character.isLowerCase(c)){
				element.append(c);
			} else if (Character.isDigit(c)){
				number.append(c);
			}
		} // end for
		
		// empty stringbuffers for the last time
		// if an element has been read
		if (element.length() != 0) {
			// lookup element
			el = et.get(element.toString());
			// if a number has been read: parse amount, enter in map and update mass.
			if (number.length() != 0) {
				amount = Integer.parseInt(number.toString());
				contents[position] = new Abundancy(element.toString(), amount);
				++position;
				mass += amount*el.getMass();
			}	else { // Insert one element
				contents[position] = new Abundancy(element.toString(), 1);
				++position;
				mass += el.getMass();
			}
		}
		Arrays.sort(contents, 0, position);
		dbe= calcDBE();
	}
	
  public double getScore() {
    return score;
  }

  public void logAndAddScore(double sc) {
    score += Math.log(sc);
    //score *= sc;
  }

  public void addScore(double sc) {
  	score += sc;
  }
  
  public void setScore(double sc) {
    score = sc;
  }

	public Vertex getVertex() {
		return vertex;
		
	}
	
	public void setVertex(Vertex vertex) {
		this.vertex = vertex;
	}

	public double getMass(){
		return mass;
	}
	
	public long getNominalMass(){
	  long nmass = 0;
	  for (Abundancy a : contents){
	    if (a != null) {
  	    try {
  	      nmass += a.amount*Math.round(elementTable.get(a.name).getMass());
  	    } catch (ElementNotFoundException e) {}
	    }
	  }
	  return nmass;
	}
	
	public double getDBE(){
		return dbe;
	}
	
	public ElementTable getElementTable(){
		return elementTable;
	}
	
	public void addElement(String toBeAdded) throws ElementNotFoundException {
		addElement(toBeAdded, 1);
	}
	
	public void addElement(String toBeAdded, int abundance) throws ElementNotFoundException {
		if (abundance == 0){
		  return;
    }
    if (abundance < 0){
		  removeElement(toBeAdded, -abundance);
      return;
    }
    mass += elementTable.get(toBeAdded).getMass()*abundance;
		int i = 0;
		while (i < position && !toBeAdded.equals(contents[i].name)) ++i;
		if (i == position) {
			contents[position] = new Abundancy(toBeAdded, abundance);
			++position;
			Arrays.sort(contents, 0, position);			
		} else {
			contents[i].amount += abundance;
		}
		dbe = calcDBE();
	}
	
	public void removeElement(String toBeRemoved) throws ElementNotFoundException {
    removeElement(toBeRemoved, 1);
	}
	
  public void removeElement(String toBeRemoved, int abundance) throws ElementNotFoundException {
    Element el = elementTable.get(toBeRemoved);
    int i = 0;
    while (i < position && !toBeRemoved.equals(contents[i].name)) ++i;
    if (i == position || contents[i].amount == 0){
      return; // or throw exception?
    } else {
      contents[i].amount -= abundance;
    }
    mass -= el.getMass()*abundance;
    Arrays.sort(contents, 0, position);
    dbe = calcDBE();
  }

  public int getElementAbundance(String which){
		int i = 0;
		while (i < position && !which.equals(contents[i].name)) ++i;
		if (i == position) return 0;
		return contents[i].amount;
	}
	
	public double getIsotopeScore() {
    return isotopeScore;
  }

  public void setIsotopeScore(double isotopeScore) {
    this.isotopeScore = isotopeScore;
  }

  public boolean isNeutral(){
    return (2*dbe % 2) == 0;
  }
  
  public boolean hasPositiveDBE(){
		/*int[] valencies = new int[5];
		for (int i = 0; i < position; ++i){
			valencies[elementTable.get(contents[i].name).getValency()] += contents[i].amount;
		}
		return ((valencies[1] - valencies[3]) % 2 == 0) &&			// filters if DBE is integer
			((valencies[3] + 2 * valencies[4] + 2) >= valencies[1]);	// filters if DBE > 0 */
		return dbe >= 0;
	}
	
	private double calcDBE() throws ElementNotFoundException {
		//int[] valencies = new int[5];
		double dbe = 1.0;
		for (int i = 0; i < position; ++i){
			//valencies[elementTable.get(contents[i].name).getValency()] += contents[i].amount;
			dbe += (0.5*elementTable.get(contents[i].name).getValency()-1.0)*contents[i].amount;
		}
		return dbe;
		//return 1.0 - valencies[1]/2.0 + valencies[3]/2.0 + valencies[4];
	}
	
	public boolean isTrueSubsetOf(Molecule parent){
        boolean wasSmaller = false;
        boolean jumpedElement = false;
        if (parent.mass <= mass) return false;
        int pIndex = 0;
        for (int i = 0; i < position; ++i){
            while (pIndex < parent.position && !contents[i].name.equals(parent.contents[pIndex].name)){
                ++pIndex;
            }
            if (pIndex>i)jumpedElement=true;
            if (pIndex >= parent.position || contents[i].amount > parent.contents[pIndex].amount) return false;
            if (contents[i].amount < parent.contents[pIndex].amount) wasSmaller = true;
        }
        // If the molecules were identical till the end of this molecule, but the parent molecule continues, return true
        if (pIndex < parent.position-1) return true;
        if (!wasSmaller) return jumpedElement;
        return wasSmaller;
	}

	public boolean isSubsetOf(Molecule parent){
		if (parent.mass < mass) return false;
		int pIndex = 0;
		for (int i = 0; i < position; ++i){
			while (pIndex < parent.position && !contents[i].name.equals(parent.contents[pIndex].name)){
				++pIndex;
			}
			if (pIndex >= parent.position || contents[i].amount > parent.contents[pIndex].amount) return false;
		}
		return true;
	}
	
	public Molecule add(Molecule summand) throws ElementNotFoundException{
		Molecule result = new Molecule(this.elementTable);
		int sIndex = 0, i = 0;
		int comparison = 0;
		while (sIndex < summand.position || i < position){
			if (sIndex >= summand.position){
				comparison = -1;
			} else if (i >= position){
				comparison = 1;
			} else {
				comparison = contents[i].compareTo(summand.contents[sIndex]);
			}
			
			if (comparison == 0){
				result.addElement(contents[i].name, contents[i].amount+summand.contents[sIndex].amount);
				++sIndex;
				++i;
			} else if (comparison < 0) {
				result.addElement(contents[i].name, contents[i].amount);				
				++i;
			} else {
				result.addElement(summand.contents[sIndex].name, summand.contents[sIndex].amount);				
				++sIndex;
			}
		}
		return result;
	}
	
	public Molecule subtract(Molecule child) throws ElementNotFoundException{
		Molecule result = new Molecule(this.elementTable);
		int cIndex = 0;
		for (int i = 0; i < position; ++i){
			while (cIndex < child.position && !contents[i].name.equals(child.contents[cIndex].name)) {
				++cIndex;
			}
			if (cIndex >= child.position){
				result.addElement(contents[i].name, contents[i].amount);
				cIndex = 0;
			}
			else result.addElement(contents[i].name, contents[i].amount-child.contents[cIndex].amount);
		}
		return result;
	}
	
	public int compareTo(Molecule another){
		return this.toString().compareTo(another.toString());
	}
	
	public boolean equals(Object o) {
		if (o == null || !(o instanceof Molecule)) return false;
		Molecule another = (Molecule) o;
		if (position != another.position) return false;
		for (int i = 0; i < position; ++i){
			if (contents[i].amount != another.contents[i].amount || !contents[i].name.equals(another.contents[i].name)) return false;
		}
		return true;
	}
	
  public int hashCode(){
    int result = 0;
    for (int i = 0; i< position; ++i){
      result ^= contents[i].amount;
    }
    return result;
  }
  
	public String noHString(){
		StringBuffer buffer = new StringBuffer();
		for (int i = 0; i < position; ++i){
			if (!contents[i].name.equals("H")){
				if (contents[i].amount > 0) buffer.append(contents[i].name);
				if (contents[i].amount > 1) buffer.append(contents[i].amount);
			}
		}
		return buffer.toString();
	}

	public String toString(){
		StringBuffer buffer = new StringBuffer();
		for (int i = 0; i < position; ++i){
			if (contents[i].amount > 0) buffer.append(contents[i].name);
			if (contents[i].amount > 1) buffer.append(contents[i].amount);
		}
		return buffer.toString();
	}
	
	public String toCmlString(){
		StringBuffer buffer = new StringBuffer();
		for (int i = 0; i < position; ++i){
			if (contents[i].amount > 0) buffer.append(contents[i].name); buffer.append(" ");
			if (contents[i].amount > 0) buffer.append(contents[i].amount); buffer.append(" ");
		}
		return buffer.toString();
	}
	
  public static Comparator<Molecule> scoreComparator(){
    return new Comparator<Molecule>() {
      public int compare(Molecule m1, Molecule m2){
        return (int) Math.signum(m2.getScore() - m1.getScore());
      }      
    };
  }
	public IsotopeDistribution calculateIsotopeDistribution(int scaled, int additional){
	  IsotopeDistribution result = new IsotopeDistribution();
	  for (Abundancy a : contents) {
	    if (a != null) {
  	    try { 
  	      result = result.fold(elementTable.get(a.name).getIsotopeDistribution().fold(a.amount));
  	    } catch (ElementNotFoundException e) {}
	    }
	  }
    int size = scaled+additional;
	  List<Peak> peaks = result.getPeaks();
	  if (peaks.size() > scaled){
	    double intSum = 0.0;
	    for (int i = 0; i<scaled; ++i){
	      intSum += peaks.get(i).getIntensity();	      
	    }
	    for (int i = 0; i<scaled; ++i){
	      peaks.get(i).setIntensity(peaks.get(i).getIntensity()/intSum);
	    }
      if (peaks.size() > size){
        for (int i = scaled; i<size; ++i){
          intSum += peaks.get(i).getIntensity();        
        }
        for (int i = scaled; i<size; ++i){
          peaks.get(i).setIntensity(peaks.get(i).getIntensity()/intSum);
        }        
      }
	    for (int i = peaks.size()-1; i>=size; --i){
	      peaks.remove(i);
	    }
	  }
	  long nomMass = this.getNominalMass();
    int i = 0;
    for (Peak p : peaks){
      p.setMass(p.getMass()+(double) nomMass +i);
      ++i;
    }
	  
	  result.setNominalMass(new Integer((int) nomMass));
	  return result;
	}
  
  public int size(){
    return position;
  }
  
  public Iterator<String> iterator(){
    return new Iterator<String>(){
      int p = 0;

      public boolean hasNext() {
        return p<position;
      }

      public String next() {
        ++p;
        return contents[p-1].name;
      }

      public void remove() {
        // TODO Auto-generated method stub
        
      }
    };
  }
}

class Abundancy implements Comparable<Abundancy>, java.io.Serializable {

  private static final long serialVersionUID = 2676221303795126054L;
  String name;
	int amount;
	
	Abundancy(String name, int amount){
		this.name = name;
		this.amount = amount;
	}
	
	public int compareTo(Abundancy theOther){
		return name.compareTo(theOther.name);
	}
	
}
