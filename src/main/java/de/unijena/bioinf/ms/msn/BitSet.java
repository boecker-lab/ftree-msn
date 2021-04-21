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
public class BitSet {
  public static long add(int element, long set){
    if (element<0 || element>63) {
      throw new IllegalArgumentException("Elements must be non-negative and <64, element given:"+element);
    }    
    return set | createMask(element);
  }
  
  public static long remove(int element, long set){
    if (element<0 || element>63) {
      throw new IllegalArgumentException("Elements must be non-negative and <64, element given:"+element);
    }    
    return set & (~createMask(element));
  }
  
  public static long singleton(int element){
    if (element<0 || element>64) {
      throw new IllegalArgumentException("Elements must be non-negative and <64, element given:"+element);
    }    
    return createMask(element);
  }
  
  public static int size(long set){
    return Long.bitCount(set);
  }
  
  public static long intersect(long set1, long set2){
    return set1 & set2; // intersection is identical to bitwise and
  }

  public static long setMinus(long minuend, long subtrahend){
    return minuend & (~subtrahend); // setMinus = andNot
  }

  public static long union(long set1, long set2){
    return set1 | set2;
  }
  
  public static boolean contains(int element, long set){
    if (element<0 || element>64) {
      return false;
    }
    return (createMask(element) & set) != 0l;
  }
  
  public static String toString(long set){
    boolean first = true;
    StringBuilder b = new StringBuilder(250);
    b.append('{');
    
    for (int i = 0; i<64; ++i){
      if ((set & 1l) == 1l){
        if (!first) {
          b.append(", ");
        } else {
          first = false;
        }
        b.append(i);
      }
      set >>>= 1;
    }
    b.append('}');
    return b.toString();    
  }
  
  private static long createMask(int pos){
    return 1l << pos;
  }
  
}
