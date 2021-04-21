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

import java.util.TreeMap;

/**
 * @author Kerstin Scheubert
 */
public class ElementTable  implements java.io.Serializable {

  private static final long serialVersionUID = -5940500125399559063L;
  private TreeMap<String, Element> contents;
	
	public ElementTable(){
		contents = new TreeMap<String, Element>();
	}
	
	public void add(Element element){
		contents.put(element.getName(), element);
	}
	
	public Element get(String name) throws ElementNotFoundException {
		Element result = contents.get(name);
		if (result == null) throw new ElementNotFoundException(name);
		return result;
	}
	
	public Element[] getAll(){
		return contents.values().toArray(new Element[contents.size()]);
	}
	
	public int size(){
		return contents.size();
	}
	
	public static ElementTable generateCHNOPS(){
		ElementTable et = new ElementTable();
		Element c = new Element("C", 12.0, 4);
		c.addIsotope(new Peak(0.0, 0.9889, 0));
		c.addIsotope(new Peak(0.003355, 0.0111, 0));
		et.add(c);
		Element h = new Element("H", 1.007825, 1);
		h.addIsotope(new Peak(0.007825, 0.99985, 0));
		h.addIsotope(new Peak(0.014102, 1.5E-4, 0));
		et.add(h);
		Element n = new Element("N", 14.003074, 3);
		n.addIsotope(new Peak(0.003074, 0.99634, 0));
		n.addIsotope(new Peak(1.09E-4, 0.00366, 0));
		et.add(n);
		Element o = new Element("O", 15.994915, 2);
		o.addIsotope(new Peak(-0.005085, 0.99762, 0));
		o.addIsotope(new Peak(-8.68E-4, 3.8E-4, 0));
		o.addIsotope(new Peak(-8.39E-4, 0.0020, 0));
		et.add(o);
		Element p = new Element("P", 30.973762, 3);
		p.addIsotope(new Peak(-0.026238, 1.0, 0));
		et.add(p);
		Element s = new Element("S", 31.972071, 2);
		s.addIsotope(new Peak(-0.027929, 0.9502, 0));
		s.addIsotope(new Peak(-0.028541, 0.0075, 0));
		s.addIsotope(new Peak(-0.032133, 0.0421, 0));
		s.addIsotope(new Peak(0.0, 0.0, 0));
		s.addIsotope(new Peak(-0.032919, 2.0E-4, 0));
		et.add(s);
		
    /*Element cl = new Element("Cl", 34.968852, 1);
    cl.addIsotope(new Peak(-0.031147, 0.7577, 0));
    cl.addIsotope(new Peak(0.0, 0.0, 0));
    cl.addIsotope(new Peak(-0.034097, 0.2423, 0));
    et.add(cl);
    
    Element f = new Element("F", 18.998403, 1);
    f.addIsotope(new Peak(-0.001597, 1.0, 0));
    et.add(f);*/
    
		return et;
	}
	
}
