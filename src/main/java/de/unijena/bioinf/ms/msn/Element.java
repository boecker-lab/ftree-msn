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

/**
 * @author Kerstin Scheubert
 */
public class Element implements Comparable<Element>, java.io.Serializable {

  private static final long serialVersionUID = 3738954610579130080L;
  public static final double ELECTRON_MASS = 0.00054858;
  public static final double PROTON_MASS = 1.00727647;
  public static final double NEUTRON_MASS = 1.00866491;
  public static final double HYDROGEN_MASS = 1.00782505;
	private String name;
	private double mass;
	private int valency;
	private IsotopeDistribution id;
	
	public Element(String name, double mass, int valency){
		if (mass < 0) System.out.println("Masses have to be positive");
		this.name = name;
		this.mass = mass;
		this.valency = valency;
		id = new IsotopeDistribution();
	}
	
	public String getName(){
	 return name;
	}

	public double getMass(){
	 return mass;
	}

	public int getValency(){
	 return valency;
	}
	
	public IsotopeDistribution getIsotopeDistribution(){
	  return id;
	}
	
	public void setName(String name){
		this.name = name;
	}

	public void setMass(double mass){
		this.mass = mass;
	}

	public void setValency(int valency){
		this.valency = valency;
	}
	
	public void setIsotopeDistribution(IsotopeDistribution id){
	  this.id = id;
	}
	
	public void addIsotope(Peak p){
	  id.addIsotope(p);
	}
	
	public String toString(){
		return name;
	}
	
	public int compareTo(Element o){
		return name.compareTo(o.name);
	}
}
