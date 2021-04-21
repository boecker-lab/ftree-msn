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

import java.util.NoSuchElementException;

/**
 * @author Kerstin Scheubert
 */
public class ColourSetsMap {
  //private FlexibleMap[] maps;
  private LongDoubleHashMap[] maps;
  private int colour;
  
  public ColourSetsMap(int numberOfColours, int colour){
    //maps = new FlexibleMap[numberOfColours-1];
    maps = new LongDoubleHashMap[numberOfColours-1];
    for (int i = 0; i<numberOfColours-1; ++i){
      //int sizeOfSet = i+2;
      //maps[i] = new FlexibleMap(numberOfColours, sizeOfSet);
      maps[i] = new LongDoubleHashMap();
    }
    this.colour = colour;
  }
  
  public void updateSet(long set, int size, double value){
    //System.out.println("Update"+BitSet.toString(set)+" "+size+" "+value);
    if (size < 2){
      throw new IllegalArgumentException("Tried to update a too small set!");
    }
    double current = maps[size-2].get(set); 
    if (Double.isNaN(current) || current < value){
      maps[size-2].put(set, value);
    }
  }

  public LongIterator keys(int size){
    //if size=1, imitate a map containing one pair with key singleton(colour) 
    if(size==1) {
      return new LongIterator() {
        private boolean called=false;
        public boolean hasNext() {
          return !called;
        }

        public long next() {
          if (called){
            throw new NoSuchElementException();
          }
          called = true;
          return BitSet.singleton(colour);
        }
        
      };
    }
    if(size==0) throw new IllegalArgumentException("Size of 0 not permitted");
    return maps[size-2].keys();
  }
  
  public DoubleIterator values(int size){
    //if size=1, imitate a map containing one pair with value 0.0 
    if(size==1) {
      return new DoubleIterator() {
        private boolean called=false;
        public boolean hasNext() {
          return !called;
        }

        public double next() {
          if (called){
            throw new NoSuchElementException();
          }
          called = true;
          return 0.0;
        }
        
      };
    }
    if(size==0) throw new IllegalArgumentException("Size of 0 not permitted");
    return maps[size-2].values();
  }
  
  public double maxScore(){
    double result = 0.0;
    for (int i = 0; i < maps.length; ++i){
      for (DoubleIterator vI = maps[i].values(); vI.hasNext();){
        double val = vI.next();
        if (val > result){
          result = val;
        }
      }
    }
    return result;
  }

  public long setWithMaxScore(){
    double score = 0.0;
    long set = BitSet.singleton(colour);
    for (int i = 0; i < maps.length; ++i){
      LongIterator kIter = maps[i].keys();
      for (DoubleIterator vIter = maps[i].values(); vIter.hasNext();){
        double val = vIter.next();
        long key = kIter.next();
        if (val > score){
          score = val;
          set = key;
        }
      }
    }
    return set;
  }

  public double get(long set, int size) {
    if (size == 1){
      if(set == BitSet.singleton(colour)){
        return 0.0;
      } else {
        return Double.NaN;
      }
    }
    if (size == 0) {
      throw new IllegalArgumentException("Cannot get value of empty set");
    }
    return maps[size-2].get(set);
  }
  
}