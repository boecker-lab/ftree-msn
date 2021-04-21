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

import java.util.ArrayList;
import java.util.List;

/**
 * @author Kerstin Scheubert
 */
public class Decomposer{
	
	private Parameters param;
	private Weights weights;
	private long[][] ERT;
	
	public Decomposer(Parameters param){
		this.param = param;
		weights = new Weights(param.et, param.precision);
		calcERT();
	}
	
  public List<Molecule> decompose(double mass, int mode){
		mass -= mode*Element.PROTON_MASS;
		//System.out.println("Decomposing: "+mass);
		double absError = param.getErrorForMass(mass);
		
		// Get MinMaxRoundingError as in decomputils.h
		double precision = weights.getPrecision(), minNegativeError = 0.0, maxPositiveError = 0.0;
		for (int i = 0; i < weights.length(); ++i) {
			
			double error = (precision * weights.getIntWeight(i)
										- weights.getWeight(i)) / weights.getWeight(i);
			if (error < 0 && error < minNegativeError) {
				minNegativeError = error;
			} else {
				if (error > 0 && error > maxPositiveError) {
					maxPositiveError = error;
				}
			}
		}
		
		// defines the range of integers to be decomposed
		long startIntegerMass = (long) Math.ceil((1 + minNegativeError) * (mass - absError) / precision);
		long endIntegerMass = (long) Math.floor((1 + maxPositiveError) * (mass + absError) / precision);
		// loops and finds decompositions for every integer mass,
		// then checks if real mass of decomposition lays in the allowed
		// error interval [mass-error; mass+error]
		List<int[]> result = new ArrayList<int[]>(), rawDecompositions;
		for (long integerMass = startIntegerMass; integerMass <= endIntegerMass; ++integerMass) {
			//System.out.println(integerMass);
			
			rawDecompositions =	integerDecompose(integerMass);
			for (int i = 0; i < rawDecompositions.size();	++i) {
				if (Math.abs(calcMass(rawDecompositions.get(i)) - mass) > absError) {
					rawDecompositions.remove(i);
					--i;
				} 
			}
			result.addAll(rawDecompositions);
			//System.out.println(toMoleculeList(result));
		}
		
		return toMoleculeList(result);
	}
	
	private void calcERT(){
		ERT = new long[(int)weights.getIntWeight(0)][weights.length()];
		long d, r, n, argmin;
		
		//Init
		ERT[0][0] = 0;
		for (int i = 1; i < ERT.length; ++i){
			ERT[i][0] = Long.MAX_VALUE; // should be infinity
		}
		
		//Filling the Table, j loops over columns
		for (int j = 1; j < ERT[0].length; ++j){
			ERT[0][j] = 0; // Init again
			d = Weights.gcd(weights.getIntWeight(0), weights.getIntWeight(j));
			for (long p = 0; p < d; p++){ // Need to start d Round Robin loops
				if (p == 0) {
					n = 0; // 0 is the min in the complete RT or the first p-loop
				} else {
					n = Long.MAX_VALUE; // should be infinity
					argmin = p;
					for (long i = p; i<ERT.length; i += d){ // Find Minimum in specific part of ERT
						if (ERT[(int)i][j-1] < n){ 
							n = ERT[(int)i][j-1];
							argmin = i;
						}
					}
					ERT[(int)argmin][j]= n;
				}
				if (n == Long.MAX_VALUE){ // Minimum of the specific part of ERT was infinity
					for (long i = p; i<ERT.length; i += d){ // Fill specific part of ERT with infinity
						ERT[(int)i][j] = Long.MAX_VALUE;
					}
				} else { // Do normal loop				
					for (long i = 1; i < ERT.length/d; ++i){ // i is just a counter
						//System.out.print("p: "+p+" ERT/d "+(ERT.length/d)+" i: "+i);
						n += weights.getIntWeight(j);
						r = n%weights.getIntWeight(0);
						if (ERT[(int)r][j-1] < n) n = ERT[(int)r][j-1]; // get the min
						ERT[(int)r][j] = n;
					}
				}
			} // end for p
		} // end for j
		
		// Output ERT for debugging
		/*for (long i = 0; i < ERT.length; ++i){
			for (int j = 0; j < ERT[0].length; ++j){
				System.out.print((ERT[(int)i][j]==Long.MAX_VALUE?"inf":Long.toString(ERT[(int)i][j])) + "\t");
			}
			System.out.println();
		}*/
	}
	
	private List<int[]> integerDecompose(long mass){
		// Find compomers
		List<int[]> result = new ArrayList<int[]>();
		int k = weights.length();
		int[] c = new int[k], deepCopy;
		long[] j = new long[k], m = new long[k], lbound = new long[k], r = new long[k];
		boolean flagWhile = false; // flag wether we are in the while-loop or not
		
		// Init
		for (int i=1; i<k; ++i){
			lbound[i] = Long.MAX_VALUE; // this is just to ensure, that lbound < m in the first iteration
		}
		
		int i = k-1;
		m[i] = mass; // m[i] corresponds to M, m[i-1] ^= m 
		while (i != k){
			if (i == 0){
				deepCopy = new int[weights.length()];
				for (int index=0; index<c.length; ++index) deepCopy[index] = c[index];
				deepCopy[0] = (int) (m[i]/weights.getIntWeight(0));
				result.add(deepCopy);
				++i; // "return" from recursion
				flagWhile = true; // in this recursion-depth we are in the while-loop, cause the next recursion (the one we just exited) was called
				m[i-1] -= weights.getLCM(i); // execute the rest of the while
				c[i] += weights.getL(i);
			} else {
				if (flagWhile){
					//System.out.println("lbound: " +lbound[i]+" m: "+m[i-1]);
					if (m[i-1] >= lbound[i]){ //currently in while loop
						//System.out.println("i: "+(i-1)+" m: "+m[i-1]);
						--i; // "do" recursive call
					} else {
						flagWhile = false; // 
					}
				} else { //we are in the for-loop
					if (j[i] < weights.getL(i) && m[i]-j[i]*weights.getIntWeight(i)>=0){
						c[i] = (int) j[i];
						m[i-1] = m[i]-j[i]*weights.getIntWeight(i);
						r[i] = m[i-1]%weights.getIntWeight(0);
						lbound[i] = ERT[(int)r[i]][i-1];
						flagWhile = true; // call the while loop
						++j[i];
					} else { //exit for loop
						// reset "function variables"
						lbound[i] = Long.MAX_VALUE;
						j[i] = 0;
						c[i] = 0;
						++i; // "return" from recursion
						if (i != k) { // only if we are not done
							flagWhile = true; // in this recursion-depth we are in the while-loop, cause the next recursion was called
							m[i-1] -= weights.getLCM(i); // execute the rest of the while
							c[i] += weights.getL(i);
						}
					}
				}
			} // end if i == 0
		} // end while
		return result;
	} // end function
	
	private double calcMass(int[] input){
		double result = 0.0;
		for (int i = 0; i < input.length; ++i){
			result += input[i]*weights.getWeight(i);
		}
		return result;
	}
	
	private List<Molecule> toMoleculeList(List<int[]> compomers){
		List<Molecule> result = new ArrayList<Molecule>(compomers.size());
		for (int i = 0; i < compomers.size(); ++i){
			try {
				result.add(new Molecule(param.et, compomers.get(i), weights));
			} 
			catch (ElementNotFoundException e) { 
				e.printStackTrace();
			}
		}
		return result;
	}
}

// Rubbish: old version of int-decomp
				/*if (j[i] < weights.getL(i) && m[i]-j[i]*a[i]>=0) { //currently in for loop, s-equal to ensure while is done the last time.
					if (m[i-1] >= lbound[i]){ //currently in while loop
						--i; // "do" recursive call
					} else { //exit while loop, execute next for loop
						if (j[i] != weights.getL(i)){ // no need to execute loop body if j == l, just do the increment
							c[i] = (int) j[i];
							m[i-1] = m[i]-j[i]*weights.getIntWeight(i);
							r2[i] = m[i-1]%weights.getIntWeight(0);
							lbound[i] = ERT[(int)r2[i]][i-1];
						}
						++j[i];
					}
				} else { // exit for loop
					// reset "function variables"
					lbound[i] = Long.MAX_VALUE;
					j[i] = 0;
					c[i] = 0;
					++i; // "return" from recursion
					// execute the rest of the while
					if (i != k) { // only if we are not done
						while-loop = true; // in this recursion-depth we are in the loops, cause the next recursion was called
						for-loop = true;
						m[i-1] -= weights.getLCM(i);
						c[i] += weights.getL(i);
					}
				}*/
