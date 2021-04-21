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

import java.util.Arrays;
import java.util.NoSuchElementException;

/**
 * Code partly taken from Javas HashMap  
 * It cannot hold 0 as key and NaN as value. 
 * Unsynchronized!
 */
/**
 * @author Kerstin Scheubert
 */
public class LongDoubleHashMap {
    
  /**
   * The default initial capacity - MUST be a power of two.
   */
  static final int DEFAULT_INITIAL_CAPACITY = 8;

  /**
   * The maximum capacity, used if a higher value is implicitly specified
   * by either of the constructors with arguments.
   * MUST be a power of two <= 1<<31.
   */
  static final int MAXIMUM_CAPACITY = 1 << 31;

  /**
   * The load factor used when none specified in constructor.
   */
  static final float DEFAULT_LOAD_FACTOR = 0.75f;

  /**
   * The table, resized as necessary. Length MUST Always be a power of two.
   */
  long[] keys;
  double[] values;

  /**
   * The number of key-value mappings contained in this map.
   */
  int size;

  /**
   * The next size value at which to resize (capacity * load factor).
   */
  int threshold;

  /**
   * The load factor for the hash table.
   */
  final float loadFactor;

  /**
   * Constructs an empty <tt>HashMap</tt> with the specified initial
   * capacity and load factor.
   *
   * @param  initialCapacity the initial capacity
   * @param  loadFactor      the load factor
   * @throws IllegalArgumentException if the initial capacity is negative
   *         or the load factor is nonpositive
   */
  public LongDoubleHashMap(int initialCapacity, float loadFactor) {
      if (initialCapacity < 0)
          throw new IllegalArgumentException("Illegal initial capacity: " +
                                             initialCapacity);
      if (initialCapacity > MAXIMUM_CAPACITY)
          initialCapacity = MAXIMUM_CAPACITY;
      if (loadFactor <= 0 || Float.isNaN(loadFactor))
          throw new IllegalArgumentException("Illegal load factor: " +
                                             loadFactor);

      // Find a power of 2 >= initialCapacity
      int capacity = 1;
      while (capacity < initialCapacity)
          capacity <<= 1;

      this.loadFactor = loadFactor;
      threshold = (int)(capacity * loadFactor);
      keys = new long[capacity];
      values = new double[capacity];
      Arrays.fill(values, Double.NaN);
 }

  /**
   * Constructs an empty <tt>HashMap</tt> with the specified initial
   * capacity and the default load factor (0.75).
   *
   * @param  initialCapacity the initial capacity.
   * @throws IllegalArgumentException if the initial capacity is negative.
   */
  public LongDoubleHashMap(int initialCapacity) {
      this(initialCapacity, DEFAULT_LOAD_FACTOR);
  }

  /**
   * Constructs an empty <tt>HashMap</tt> with the default initial capacity
   * (8) and the default load factor (0.75).
   */
  public LongDoubleHashMap() {
      this.loadFactor = DEFAULT_LOAD_FACTOR;
      threshold = (int)(DEFAULT_INITIAL_CAPACITY * DEFAULT_LOAD_FACTOR);
      keys = new long[DEFAULT_INITIAL_CAPACITY];
      values = new double[DEFAULT_INITIAL_CAPACITY];
      for (int i = 0; i<values.length; ++i){
        values[i] = Double.NaN;
      }
  }

  // internal utilities

  /**
   * Applies a supplemental hash function to a given hashCode, which
   * defends against poor quality hash functions.  This is critical
   * because HashMap uses power-of-two length hash tables, that
   * otherwise encounter collisions for hashCodes that do not differ
   * in lower bits. Note: Null keys always map to hash 0, thus index 0.
   */
  static int hash(long key, int tries) {
      // This function ensures that hashCodes that differ only by
      // constant multiples at each bit position have a bounded
      // number of collisions (approximately 8 at default load factor).
      int h = (int) (key^(key>>>32));
      //h ^= (h >>> 20) ^ (h >>> 12);
      //h ^= (h >>> 7) ^ (h >>> 4);
      return (int) (h+0.5*tries+0.5*tries*tries);
  }

  /**
   * Returns index for hash code h.
   */
  static int indexFor(int h, int length) {
      return h & (length-1);
  }

  /**
   * Returns the number of key-value mappings in this map.
   *
   * @return the number of key-value mappings in this map
   */
  public int size() {
      return size;
  }

  /**
   * Returns <tt>true</tt> if this map contains no key-value mappings.
   *
   * @return <tt>true</tt> if this map contains no key-value mappings
   */
  public boolean isEmpty() {
      return size == 0;
  }

  /**
   * Returns the value to which the specified key is mapped,
   * or {@code null} if this map contains no mapping for the key.
   *
   * <p>More formally, if this map contains a mapping from a key
   * {@code k} to a value {@code v} such that {@code (key==null ? k==null :
   * key.equals(k))}, then this method returns {@code v}; otherwise
   * it returns {@code null}.  (There can be at most one such mapping.)
   *
   * <p>A return value of {@code null} does not <i>necessarily</i>
   * indicate that the map contains no mapping for the key; it's also
   * possible that the map explicitly maps the key to {@code null}.
   * The containsKey operation may be used to
   * distinguish these two cases.
   *
   *  see put(Object, Object)
   */
  public double get(long key) {
      if (key == 0){
        throw new IllegalArgumentException("0 is never a key of this map");
      }
      int index = indexFor(hash(key, 0), keys.length);;
      long k = keys[index];
      for (int tries = 1; k != 0 && k != key && tries < keys.length; ++tries){
        index = indexFor(hash(key, tries), keys.length);
        k = keys[index];
      }
      if (k==key){
        return values[index];
      }
      return Double.NaN;
  }
  /**
   * Associates the specified value with the specified key in this map.
   * If the map previously contained a mapping for the key, the old
   * value is replaced.
   *
   * @param key key with which the specified value is to be associated
   * @param value value to be associated with the specified key
   * @return the previous value associated with <tt>key</tt>, or
   *         <tt>null</tt> if there was no mapping for <tt>key</tt>.
   *         (A <tt>null</tt> return can also indicate that the map
   *         previously associated <tt>null</tt> with <tt>key</tt>.)
   */
  public void put(long key, double value) {
    if (key == 0){
      throw new IllegalArgumentException("0 may not be entered into this map");
    }
    if (Double.isNaN(value)){
      throw new IllegalArgumentException("NaN may not be entered into this map");
    }
    
    int index = indexFor(hash(key, 0), keys.length);;
    long k = keys[index];
    for (int tries = 1; k != 0 && k != key && tries < keys.length; ++tries){
      index = indexFor(hash(key, tries), keys.length);
      k = keys[index];
    }
    if (k != 0 && k!= key){
      throw new RuntimeException("Map is full");
    }
    keys[index] = key;
    values[index] = value;
    ++size;
    if (size >= threshold){
      resize(2 * keys.length);
    }
  }

  /**
   * Rehashes the contents of this map into a new array with a
   * larger capacity.  This method is called automatically when the
   * number of keys in this map reaches its threshold.
   *
   * If current capacity is MAXIMUM_CAPACITY, this method does not
   * resize the map, but sets threshold to Integer.MAX_VALUE.
   * This has the effect of preventing future calls.
   *
   * @param newCapacity the new capacity, MUST be a power of two;
   *        must be greater than current capacity unless current
   *        capacity is MAXIMUM_CAPACITY (in which case value
   *        is irrelevant).
   */
  void resize(int newCapacity) {
      long[] oldKeys = keys;
      double[] oldValues = values;
      int oldCapacity = oldKeys.length;
      if (oldCapacity == MAXIMUM_CAPACITY) {
          threshold = Integer.MAX_VALUE;
          return;
      }

      keys = new long[newCapacity];
      values = new double[newCapacity];
      Arrays.fill(values, Double.NaN);
      size = 0; //size will be set by put
      threshold = (int)(newCapacity * loadFactor);
      for (int i = 0; i<oldKeys.length; ++i){
        if (oldKeys[i] != 0){
          put(oldKeys[i], oldValues[i]);
        }
      }
  }

    public static class KeyIterator implements LongIterator{
      private int index;
      private LongDoubleHashMap map;
      
      public KeyIterator(LongDoubleHashMap map){
        this.map = map;
        index=0;
        if (map.size > 0){
          while (index < map.keys.length && map.keys[index] == 0){
            ++index;
          }
        } else {
          index = map.keys.length;
        }
      }
      
      public long next(){
        if (index >= map.keys.length){
          throw new NoSuchElementException();
        }
        long result = map.keys[index];
        ++index;
        while (index < map.keys.length && map.keys[index] == 0){
          ++index;
        }
        return result;
      }
      
      public boolean hasNext(){
        return index < map.keys.length;
      }
      
    }

    public static class ValueIterator implements DoubleIterator{
      private int index;
      private LongDoubleHashMap map;
      
      public ValueIterator(LongDoubleHashMap map){
        this.map = map;
        index=0;
        if (map.size > 0){
          while (index < map.values.length && Double.isNaN(map.values[index])){
            ++index;
          }
        } else {
          index = map.values.length;
        }
      }
      
      public double next(){
        if (index >= map.values.length){
          throw new NoSuchElementException();
        }
        double result = map.values[index];
        ++index;
        while (index < map.values.length && Double.isNaN(map.values[index])){
          ++index;
        }
        return result;
      }
      
      public boolean hasNext(){
        return index < map.values.length;
      }
      
    }

  public LongIterator keys() {
    return new KeyIterator(this);
  }

  public DoubleIterator values() {
    return new ValueIterator(this);
  }

}
