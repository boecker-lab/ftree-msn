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

import java.lang.Math;
import java.util.Arrays;

/**
 * @author Kerstin Scheubert
 */
public class Weights {

  private Weight[] weights;
  private double precision;
  
  public Weights(ElementTable et, double precision){
  	this.precision = precision;
  	Element[] elements = et.getAll();
  	weights = new Weight[elements.length];
  	for (int i = 0; i < elements.length; ++i){
  		weights[i] = new Weight();
  		weights[i].name = elements[i].getName();
  		weights[i].weight = elements[i].getMass();
  	}
  	InitIntWeights();
  }
  
  /*public Weights(String[] names, double[] weights, double precision){
  	// better: throw exception
    if (names.length != weights.length) System.out.println("Error: names and weights do not match");
    this.names = names;
    this.weights = weights;
    this.precision = precision;
    InitIntWeights();
  }*/
  
  public double getPrecision(){
  	return precision;
  }
  
  public String getName(int index){
  	return weights[index].name;
  }
  
  public double getWeight(int index){
  	return weights[index].weight;
  }
  
  public double getWeight(String name){
  	for (int i = 0; i<weights.length; ++i){
  		if (weights[i].name == name) return weights[i].weight;
  	}
  	// better: throw exception
    System.out.println("Error: Name not found");  	
  	return 0.0;
  }

  public long getIntWeight(int index){
  	return weights[index].intWeight;
  }
  
  public long getIntWeight(String name){
  	for (int i = 0; i<weights.length; ++i){
  		if (weights[i].name == name) return weights[i].intWeight;
  	}
  	// better: throw exception
    System.out.println("Error: Name not found");  	
  	return 0;
  }
  
  public long getL(int i){
  	return weights[i].l;
  }
  
  public long getLCM(int i){
  	return weights[i].lcm;
  }

  public int length(){
  	return weights.length;
  }

	public static long gcd(long u, long v) {
		long r = 0;

		while (v != 0) {
			r = u % v;
			u = v;
			v = r;
		}
		return u;
	}
	
  private void InitIntWeights(){
		for (int i = 0; i < weights.length; ++i) {
			weights[i].intWeight = Math.round(weights[i].weight / precision);
		}
		divideByGCD();
		Arrays.sort(weights);
		calcLCMs();
  }
  
  private void divideByGCD(){
    if (weights.length <= 1) return;
		long d = gcd(weights[0].intWeight, weights[1].intWeight);
		for (int i = 2; i < weights.length; ++i) {
			d = gcd(d, weights[i].intWeight);
			if (d == 1) {
				return;
			}
		}
		// if we're here: d != 1
		precision *= d;
	
		// rescales the integer weights.
		for (int i = 0; i < weights.length; ++i) {
			weights[i].intWeight /= d;
		}
  }
  
  private void calcLCMs(){
  	weights[0].l = 1;
  	weights[0].lcm = weights [0].intWeight;
  	for (int i = 1; i < weights.length; ++i){
  		// lcm_i = a_0 * a_i / gcd(a_0, a_i)
  		// l_i = lcm_i / a_i
  		weights[i].l = weights[0].intWeight/gcd(weights[0].intWeight, weights[i].intWeight);
  		weights[i].lcm = weights[i].l*weights[i].intWeight;
  	}
  }
}

class Weight implements Comparable{
  String name;
  double weight;
  long intWeight;
	long l;
  long lcm;
	
	public int compareTo(Object anotherWeight){
		Weight w2 = (Weight) anotherWeight;
		if (this.weight < w2.weight) return -1;
		if (this.weight == w2.weight) return 0;
		return 1;
	}
}
