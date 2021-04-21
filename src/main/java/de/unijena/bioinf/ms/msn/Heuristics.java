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

import java.util.List;
import java.util.LinkedList;
import java.util.Collections;

/**
 * @author Kerstin Scheubert
 */
public class Heuristics{
	public static void MST(Vertex PMD, List<Edge> edges, int peakCount, int maxDecomp){
    Collections.sort(edges);
    LinkedList<Edge> result = new LinkedList<Edge>();
    double score = 0.0;
    
		boolean[][] used = new boolean[peakCount][maxDecomp]; // activeVertices.length = #colours
    
		UnionFind<Molecule> uf = new UnionFind<Molecule>();
		for (Edge e : edges){
			if (!used[e.toColour()][e.toNumber()] && uf.find(e.from().getLabel()) != uf.find(e.to().getLabel())){
					uf.union(e.from().getLabel(), e.to().getLabel());
					used[e.toColour()][e.toNumber()] = true;
					result.add(e);
					score += e.getScore();
			}
    }
    
  	PMD.addScore(score);
  	PMD.setOptimalTree(result);    
	}
	
	public static void greedy(Vertex PMD, List<Edge> edges, int peakCount, int numberCount){
		int[] colours = new int[peakCount];
    boolean[][] used = new boolean[peakCount][numberCount];		
    Collections.sort(edges);
    LinkedList<Edge> result = new LinkedList<Edge>();
    double score = 0.0;
		
		for (int i = 0; i < colours.length; ++i){
      colours[i] = -1;
    }

		for (Edge e : edges){
			boolean toOK = (colours[e.toColour()] == -1 || colours[e.toColour()] == e.toNumber()),
							fromOK = (colours[e.fromColour()] == -1 || colours[e.fromColour()] == e.fromNumber());
			if (!used[e.toColour()][e.toNumber()] && toOK && fromOK){
				colours[e.toColour()] = e.toNumber();
				colours[e.fromColour()] = e.fromNumber();
				used[e.toColour()][e.toNumber()] = true;
				result.add(e);
				score += e.getScore();			
			}
		}
  	PMD.addScore(score);
  	PMD.setOptimalTree(result);    		
	}
	
	public static void topDown(Vertex PMD, int peakCount, int numberCount){
		Vertex current = PMD;
		int[] colours = new int[peakCount];
    boolean[][] used = new boolean[peakCount][numberCount];		
    LinkedList<Edge> result = new LinkedList<Edge>();
    double score = 0.0;
    boolean done = false, change = false;
    
    while (!done){
    	Collections.sort(current.outEdges());
    	change = false;
    	for (Edge e : current.outEdges()){
				boolean toOK = (colours[e.toColour()] == -1 || colours[e.toColour()] == e.toNumber()),
								fromOK = (colours[e.fromColour()] == -1 || colours[e.fromColour()] == e.fromNumber());
				if (!used[e.toColour()][e.toNumber()] && toOK && fromOK){
					colours[e.toColour()] = e.toNumber();
					colours[e.fromColour()] = e.fromNumber();
					used[e.toColour()][e.toNumber()] = true;
					result.add(e);
					score += e.getScore();
					current = e.to();
					change = true;
					break;
				}    		
    	}
    	if (!change) {
    		if (current == PMD) done = true;
    		else current = PMD;
    	}
    }
  	PMD.addScore(score);
  	PMD.setOptimalTree(result);    		    
    
	}
}
