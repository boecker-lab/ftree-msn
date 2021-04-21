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
public class FlexibleMap {
  private LongDoubleHashMap hashMap;
  private double[] array;
  private static int[] binom;
  private final int setSize;
  private final int largestSetSize;
  
  // Testing
  public static void main(String[] args){
    FlexibleMap map = new FlexibleMap(3, 12);
    for (long i= 1; i<=1000; ++i){
      if (Long.bitCount(i) == 3){
        map.put(i, i);
      }
    }
    System.out.println("put done");
    for (long i= 1; i<=100; ++i){
      double val = map.get(i);
      if (!Double.isNaN(val) && val != i){
        System.out.println(val+"!="+i);
        //return;
      }
    }    
    /*for (long k : map.keys){
      System.out.print(k+" ");
    }
    System.out.println();

    for (double v : map.values){
      System.out.print(v+" ");
    }
    System.out.println();*/
    
    /*LongIterator kIter = map.keys();
    for(DoubleIterator vIter = map.values(); vIter.hasNext();){
      System.out.println(kIter.next()+" "+vIter.next());
    }
    if (kIter.hasNext()){
      System.out.println("Still keys...");
    }*/
  }

  FlexibleMap(int setSize, int largestSetSize){
    if (binom == null){
      init(largestSetSize);
    }
    if (setSize > 2){
      hashMap = new LongDoubleHashMap();
    } else {
      array = new double[binom(largestSetSize, setSize)];
      for (int i = 0; i < array.length; i++) {
        array[i] = Double.NaN;
      }
    }
    this.setSize = setSize;
    this.largestSetSize = largestSetSize;
  }
  
  public double get(long key){
    if(hashMap != null){
      return hashMap.get(key);
    }
    return array[getArrayIndex(key)];
  } 
  
  public void put(long key, double value){
    if(hashMap != null){
      hashMap.put(key, value);
      checkSwitch();
      return;
    }
    array[getArrayIndex(key)]=value;
  }
  
  public LongIterator keys(){
    if(hashMap != null){
      return hashMap.keys();
    }
    return new KeyIterator(this);
  }
  
  public static class KeyIterator implements LongIterator{
    private int index;
    private FlexibleMap map;
    
    public KeyIterator(FlexibleMap map){
      index=0;
      this.map = map;
      while (index < map.array.length && Double.isNaN(map.array[index])){
        ++index;
      }
    }

    public boolean hasNext() {
      return index < map.array.length;
    }

    public long next() {
      if (index >= map.array.length){
        throw new NoSuchElementException();
      }
      long result = map.getSetFromIndex(index);
      ++index;
      while (index < map.array.length && Double.isNaN(map.array[index])){
        ++index;
      }      
      return result;
    }    
  }
  

  public DoubleIterator values(){
    if(hashMap != null){
      return hashMap.values();
    }
    return new ValueIterator(this);
  }
  
  public static class ValueIterator implements DoubleIterator{
    private int index;
    private FlexibleMap map;
    
    public ValueIterator(FlexibleMap map){
      index=0;
      this.map = map;
      ++index;
      while (index < map.array.length && Double.isNaN(map.array[index])){
        ++index;
      }
    }

    public boolean hasNext() {
      return index < map.array.length;
    }

    public double next() {
      if (index >= map.array.length){
        throw new NoSuchElementException();
      }
      double result = map.array[index];
      while (index < map.array.length && Double.isNaN(map.array[index])){
        ++index;
      }      
      return result;
    }    
  }

  private int getArrayIndex(long set){
    int arrayPos = 0, elementsSeen=0;
    for (int setpos = 0; setpos<largestSetSize && elementsSeen<setSize; ++setpos){
      if ((set & 1l)== 1l){
        ++elementsSeen;
      } else {
        arrayPos += binom(largestSetSize-setpos-1,setSize-elementsSeen-1);
      }
      set >>>= 1;
    }
    return arrayPos;
  }
  
  private long getSetFromIndex(int index){
    long set = 0;
    int elementsWritten = 0, checkval = 0, setpos = 0;
    for (; setpos<largestSetSize && elementsWritten<setSize; ++setpos){
      checkval = binom(largestSetSize-setpos-1,setSize-elementsWritten-1);
      //checks whether set falls in the current range of the index...
      if (index <= checkval){
        set &= 1l;
        ++elementsWritten;
        // if index == checkval fill the remaining ones
        if (index == checkval){
          while (elementsWritten<setSize){
            set <<= 1;
            set &= 1l;
            ++elementsWritten;            
          }
        }
      } else {
        index -= checkval;
      }
      set <<= 1;
    }
    //if all elements written, shift a suitable number of zeros into the set.
    if (setpos < largestSetSize){
      set <<= (largestSetSize-setpos);
    }
    return set;
  }
  
  private void checkSwitch(){
    if(hashMap == null){
      return;
    }
    if (hashMap.size() >= binom(largestSetSize, setSize)/2){
      switchToArray();
    }
  }
  
  private void switchToArray(){
    System.out.println("Array");
    array = new double[binom(largestSetSize, setSize)];
    for (int i = 0; i < array.length; i++) {
      array[i] = Double.NaN;
    }    
    LongIterator k = hashMap.keys();
    for (DoubleIterator v = hashMap.values(); v.hasNext() || k.hasNext();){
      array[getArrayIndex(k.next())] = v.next();
    }
    hashMap = null;
  }
  
  //Binomial coefficient stuff:
  private static void init(int n){
    int nHalbe = (int) Math.floor(n/2.0),
        corr = n%2==1?0:1;
    binom = new int[nHalbe*(nHalbe+1-corr)]; //(summe 1+2 ... n/2)*2-1 (-n/2 falls n gerade)
  }
  
  private static int binom(int n, int k){
    if (k > n/2.0){
      k = n-k;
    }
    if (k < 0){
      return 0;
    }
    if (k== 0){
      return 1;
    }
    int nHalbe = (int) Math.floor((n-1)/2.0),
        index = nHalbe*(nHalbe+1-(n%2)) + k;
    //System.out.println("binom"+n+" "+k+" "+index);
    if (binom[index] != 0){
      return binom[index];
    }
    int result = 1;
    for (int i=1; i<k; ++i){
      result *= n+1-i;
      result /= i;
    }
    binom[index]= result;
    return result;
  }
}