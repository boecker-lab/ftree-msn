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

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * @author Kerstin Scheubert
 */
public class Tree extends Graph {
  private Vertex orgRoot;
  private Vertex root;
  
  public Tree(List<Edge> edges){
    super(new LinkedList<Vertex>(), new LinkedList<Edge>());
    root = null;
    orgRoot = null;
    Map<Vertex, Vertex> newVertices = new HashMap<Vertex, Vertex>(edges.size()+1);
    Vertex newFrom, newTo, from , to;
    Edge newEdge;
    for(Edge e : edges){
      from = e.from();
      newFrom = newVertices.get(from);
      if (newFrom == null){
        newFrom = new Vertex(from.getLabel(), from.getPeak(), from.getColour(), from.getNumber(), from.isPMD());
        newVertices.put(from, newFrom);
        getVertices().add(newFrom);
      }
      to = e.to();
      newTo = newVertices.get(to);
      if (newTo == null){
        newTo = new Vertex(to.getLabel(), to.getPeak(), to.getColour(), to.getNumber(), to.isPMD());
        newVertices.put(to, newTo);
        getVertices().add(newTo);
      }
      newEdge = new Edge(newFrom, newTo, e.getScore());
      getEdges().add(newEdge);
      newFrom.addOutEdge(newEdge);
      newTo.addInEdge(newEdge);      
    }
    
    // search and set Roots, check for Tree property!
    for(Vertex v : getVertices()){
      if (v.inEdges().size() > 1){
        throw new RuntimeException("Tree property violated");
      }
      if (v.inEdges().size() == 0){
        if (root == null){
          root = v;
          for (Map.Entry<Vertex, Vertex> entry : newVertices.entrySet()){
            if (entry.getValue() == root){
              orgRoot = entry.getKey();
            }
          }
        } else {
          throw new RuntimeException("Tree has more than one connected component");
        }
      }
    }
  }

  public Vertex getOrgRoot() {
    return orgRoot;
  }

  public Vertex getRoot() {
    return root;
  }
  
  public void addEdge(Edge e){
    if (e.to() == root){
      root = e.from();
      orgRoot = null;
      getEdges().add(e);
      getVertices().add(e.from());
      e.to().addInEdge(e);
      e.from().addOutEdge(e);
    } else {
      if (getVertices().contains(e.to())){
        throw new RuntimeException("Treeproperty violated at Vertex "+e.to());
      }
      getEdges().add(e);
      getVertices().add(e.to());
      e.to().addInEdge(e);
      e.from().addOutEdge(e);      
    }
  }

  public void addEdgeRemoveCycle(Edge e) {
    if (e.to().inEdges().isEmpty()){
      addEdge(e);
    } else {
      Edge concurrEdge = e.to().inEdges().get(0);
      getEdges().remove(concurrEdge);
      //TODO: the following cures the asparagine bug
      concurrEdge.from().outEdges().remove(concurrEdge);
      e.to().inEdges().clear();
      getEdges().add(e);
      e.to().addInEdge(e);
      e.from().addOutEdge(e);
    }
  } 
  
	public double scoreColoursUpstairs(Vertex pv,Vertex cv){
		double sum=pv.getPeak().scoreColours(cv.getPeak());
		if(pv!=root)sum+=scoreColoursUpstairs(pv.inEdges().get(0).from(),cv);
		return sum;
	}
	
	public double scoreColoursDownstairs(Vertex pv,Vertex cv){
		double sum=pv.getPeak().scoreColours(cv.getPeak());
		for(Edge e:cv.outEdges()){
			sum+=scoreColoursDownstairs(pv,e.to());
		}
		return sum;
	}
	
	public double scoreColours(){
		return scoreColours(root);			
	}

	private double scoreColours(Vertex root){
	 	double sum=scoreColoursDownstairs(root,root);
		for(Edge e:root.outEdges()){
			sum+=scoreColours(e.to());
		}
		return sum;
	}
	
	public String dotStringTransitiveClosure() {
		return dotStringTransitiveClosure(root);
	}
	
	private StringBuffer dotStringDownstairs(Vertex pv,Vertex cv){
		if(pv==null||cv==null)return new StringBuffer("");
		StringBuffer result=new StringBuffer();
		double value=pv.getPeak().scoreColours(cv.getPeak());		
		if(value!=0){
			result.append(pv.getLabel()+" -> " + cv.getLabel()+"[style=dotted,label=\"Score: "+value+"\"]\n");
		}
		for(Edge e:cv.outEdges()){
			result.append(dotStringDownstairs(pv,e.to()));
		}
		return result;
	}
	
	private String dotStringTransitiveClosure(Vertex v) {
		if(v==null)return "";
		StringBuffer result=new StringBuffer(dotStringDownstairs(v,v));
		for(Edge e:v.outEdges()){
			result.append(dotStringTransitiveClosure(e.to()));
		}
		return result.toString();
	}
	
	public String dotStringMS3Score() {
		StringBuffer result=new StringBuffer();
		for(Vertex pv:getVertices()){
			for(Vertex cv:getVertices()){
				double value=pv.getPeak().scoreColours(cv.getPeak());		
				if(value!=0){
					result.append(pv.getLabel()+" -> " + cv.getLabel()+"[style=dotted,label=\"Score: "+value+"\"]\n");
				}

			}
		}
		return result.toString();
	}
}
