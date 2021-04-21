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

import gurobi.GRBEnv;
import gurobi.GRBException;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.LinkedList;
import java.util.Collections;
import java.util.Map;
import java.util.Map.Entry;

import net.sf.javailp.Constraint;
import net.sf.javailp.Linear;
import net.sf.javailp.OptType;
import net.sf.javailp.Problem;
import net.sf.javailp.Result;
import net.sf.javailp.ResultImpl;
import net.sf.javailp.Solver;
import net.sf.javailp.SolverFactory;
import net.sf.javailp.SolverFactoryGurobi;


/**
 * @author Kerstin Scheubert
 */
public class Graph{

	private List<Vertex> v;
	private List<Edge> e;

	public Graph(List<Vertex> vertices, List<Edge> edges){
		v = vertices;
		e = edges;
	}

	public List<Vertex> getVertices(){
		return v;
	}

	public List<Edge> getEdges(){
		return e;
	}

	public static boolean ilpsolver(Vertex PMD, Graph g, int colors, Parameters param, boolean noTransitiveClosureScore){	//naiv	
				SolverFactory factory = new SolverFactoryGurobi();
				System.out.print(PMD+" initalizing... ");
				factory.setParameter(Solver.VERBOSE, 0);
				//	        factory.setParameter(Solver.TIMEOUT, 100); // set timeout to 100 seconds

				Problem problem = new Problem();

				//transitive closure        
				HashMap<Integer,HashMap<Integer,Double>> ms3score= PMD.getPeak().ms3score;
				HashMap<Integer,HashMap<Integer,Double>> colorMap=new HashMap<Integer,HashMap<Integer,Double>>();
				for(Edge e:g.getEdges()){
					double score=0;
					if(ms3score.containsKey(e.from().getColour())&&ms3score.get(e.from().getColour()).containsKey(e.to().getColour())){
						score=ms3score.get(e.from().getColour()).get(e.to().getColour());	        		
					}	        	
					if(!colorMap.containsKey(e.from().getColour())){
						colorMap.put(e.from().getColour(), new HashMap<Integer,Double>());
					}	        
					colorMap.get(e.from().getColour()).put(e.to().getColour(), score);	        	
				}

				List<TransitiveEdgeColor> transitiveEdges=new ArrayList<TransitiveEdgeColor>();
				List<EdgeColor> edgesColors=new ArrayList<EdgeColor>();
				for(Entry<Integer,HashMap<Integer,Double>> entry1:colorMap.entrySet()) {
					for(Entry<Integer,Double> entry2:entry1.getValue().entrySet()) {
						transitiveEdges.add(new TransitiveEdgeColor(entry1.getKey(),entry2.getKey(),entry2.getValue()));
						edgesColors.add(new EdgeColor(entry1.getKey(),entry2.getKey()));
					}	        	
				}	       

				//OBJ
				Linear linear = new Linear();
				for( Edge e : g.getEdges()) {	        	
					linear.add(e.getScore(), e);
				}	  
				//transitive closure
				for(TransitiveEdgeColor t : transitiveEdges) {	        	
					linear.add(t.score, t);
				}

				problem.setObjective(linear, OptType.MAX);

				// Constraints
				for(Vertex v : g.getVertices()) {
					Linear l1 = new Linear();
					for( Edge uv : v.inEdges()) {
						l1.add(1, uv);
					}
					problem.add(l1, "<=", 1);  //(1)	         

					if(v!=PMD) {
						for(Edge vw : v.outEdges()) {
							Linear l2 = new Linear(l1);
							l2.add(-1,vw);
							problem.add(l2, ">=", 0); //(2)	       
						}
					}
				}
				for(int i=0; i< colors; i++) {
					Linear l3 = new Linear();
					for(Vertex v : g.getVertices()) {
						if(v.getColour()==i) {
							for(Edge uv : v.inEdges()) {
								l3.add(1,uv);
							}
						}
					}
					problem.add(l3, "<=", 1); //(3)	        
				}	 

				//transitive closure
				for(Edge e:g.getEdges()){
					for(EdgeColor ec:edgesColors){
						if(e.from().getColour()==ec.from()&&e.to().getColour()==ec.to()){
							Linear l4=new Linear();
							l4.add(1,e);
							l4.add(-1,ec);
							problem.add(l4,"<=",0);	        			
						}
					}	        		
				}

				for(EdgeColor ec:edgesColors){
					Linear l5=new Linear();
					for(Edge e:g.getEdges()){
						if(e.from().getColour()==ec.from()&&e.to().getColour()==ec.to()){	        		
							l5.add(-1,e);	        				        
						}
					}	 
					l5.add(1,ec);
					problem.add(l5,"<=",0);
				}	        	       

				for(EdgeColor ebc:edgesColors){	
					for(TransitiveEdgeColor tac:transitiveEdges){	        	
						if(tac.to()==ebc.to()){
							boolean found=false;
							for(TransitiveEdgeColor tab:transitiveEdges){	        		       	        	
								if(tac.from()==tab.from()&&tab.to()==ebc.from()){
									found=true;
									Linear l6=new Linear();
									l6.add(-1,tac);
									l6.add(1,tab);
									l6.add(1,ebc);
									problem.add(l6,"<=",1);

									Linear l7=new Linear();
									l7.add(1,tac);
									l7.add(1,ebc);
									l7.add(-1,tab);
									problem.add(l7,"<=",1);
								}
							}	        		
							if(!found&&tac.from()!=ebc.from()){
								Linear l6=new Linear();
								l6.add(-1,tac);
								l6.add(1,ebc);
								problem.add(l6,"<=",1);

								Linear l7=new Linear();
								l7.add(-1,tac);
								l7.add(-1,ebc);
								problem.add(l7,">=",-1);
							}
						}	        		
					}
				}

				for(EdgeColor eab:edgesColors){
					for(TransitiveEdgeColor tab:transitiveEdges){
						if(eab.from()==tab.from()&&eab.to()==tab.to()){
							Linear l8=new Linear();
							l8.add(1,eab);
							l8.add(-1,tab);
							problem.add(l8,"<=",0);
						}
					}
				}

				for(TransitiveEdgeColor tab:transitiveEdges){
					Linear l9=new Linear();
					for(EdgeColor ecb:edgesColors){
						if(tab.to()==ecb.to()){
							l9.add(-1,ecb);
						}
					}
					l9.add(1,tab);
					problem.add(l9,"<=",0);
				}

				for(Edge e : g.getEdges()) {
					problem.setVarType(e, Boolean.class);
				}

				//transitive closure
				for(TransitiveEdgeColor t : transitiveEdges) {
					problem.setVarType(t, Boolean.class);
				}
				for(EdgeColor e:edgesColors){
					problem.setVarType(e, Boolean.class);
				}

				Solver solver = factory.get(); // you should use this solver only once for one problem
				System.out.print("starting... ");
				Result ilpresult = solver.solve(problem);			
				
				if(ilpresult!=null){
					System.out.print(" solved.\n");
					LinkedList<Edge> tree = new LinkedList<Edge>();
					for(Edge e : g.getEdges()) {
						if(ilpresult.get(e).intValue() > 0) {
							tree.add(e);
						}
					}
			
					PMD.addScore(ilpresult.getObjective().doubleValue());	        
					PMD.setOptimalTree(tree);
					PMD.optimalTree(noTransitiveClosureScore);
					return true;
				}else{
					System.out.print(" failed.\n");
					return false;
				}
	}
	
	public static void scoreGraphMST(Vertex PMD, Graph g, int[] amountOfColour, int[][] mapping, int maxDecomp, boolean noTransitiveClosureScore){

		long callnumber = 1;
		for (int a : amountOfColour){
			if (a != 0) {
				callnumber *= a;
			}
		}
		//System.out.println("Calls: "+callnumber);

		int[] colournumbers = new int[amountOfColour.length];
		int cI = 0; //level on which we currently change the colour
		AlgoResult result, bestResult = null;

		while (cI != -1){
			if (colournumbers[cI] < amountOfColour[cI] || (colournumbers[cI] == 0 && amountOfColour[cI] == 0)){
				if (cI == colournumbers.length-1) {
					result = findMaxTree(g, colournumbers, mapping);
					//result = mst(edgesList, colournumbers, mapping, maxDecomp);
					if (bestResult == null || result.score > bestResult.score){
						bestResult = result;
					}
					++colournumbers[cI];
				} else {
					++cI;
				}
			} else {
				colournumbers[cI] = 0;
				--cI;
				if (cI >= 0) ++colournumbers[cI];
			}
		}

		PMD.addScore(bestResult.score);
		//System.out.println(PMD+" "+bestResult.graph);
		PMD.setOptimalTree(bestResult.graph);
		PMD.optimalTree(noTransitiveClosureScore);
	}

	private static AlgoResult findMaxTree(Graph g, int[] activeVertices, int[][] mapping) {
		List<Vertex> vertices = g.getVertices();
		Collections.sort(vertices); // Sort vertices inverse topologically (bottom-up)
		Map<Edge, Double> scores = new HashMap<Edge, Double>(g.getEdges().size());
		List<Edge> result = new LinkedList<Edge>();
		double resultScore = 0.0;

		for (Edge e : g.getEdges()){
			if (isInCurrGraph(e, activeVertices, mapping)){
				scores.put(e, e.getScore());
			}
		}
		for (Vertex v : vertices){
			double sum = 0.0;
			for (Edge e : v.outEdges()){
				if (isInCurrGraph(e, activeVertices, mapping)){
					sum += e.getScore();
				}
			}
			for (Edge e : v.inEdges()){
				if (isInCurrGraph(e, activeVertices, mapping)){
					scores.put(e, e.getScore()+sum);
				}
			}      
		}

		Collections.reverse(vertices); //Now move top-down
		Map<Vertex, Boolean> reachable = new HashMap<Vertex, Boolean>(vertices.size());
		for (Vertex v : vertices){
			reachable.put(v, true);
		}

		for (Vertex v : vertices){
			double max = Double.NEGATIVE_INFINITY, score = Double.NEGATIVE_INFINITY;
			Edge maxEdge = null;
			for (Edge e : v.inEdges()){
				if (isInCurrGraph(e, activeVertices, mapping)){
					score = scores.get(e);
					if (score > max){
						max = score;
						maxEdge = e;
					}
				}
			}
			if (max > 0){
				result.add(maxEdge);
				resultScore += maxEdge.getScore();
			}
		}
		return new AlgoResult(resultScore, result);
	}

	private static boolean isInCurrGraph(Edge e, int[] activeVertices, int[][] mapping){
		return e.fromNumber() == mapping[e.fromColour()][activeVertices[e.fromColour()]] &&
		e.toNumber() == mapping[e.toColour()][activeVertices[e.toColour()]];
	}

	/*private static AlgoResult mst(List<Edge> edges, int[] activeVertices, int[][] mapping, int maxDecomp){
    Collections.sort(edges);
    List<Edge> result = new LinkedList<Edge>();
    double score = 0.0;

		boolean[][] used = new boolean[activeVertices.length][maxDecomp]; // activeVertices.length = #colours

		UnionFind<Vertex> uf = new UnionFind<Vertex>();
		for (Edge e : edges){
      if (e.getScore() < 0) break;
			if (isInCurrGraph(e, activeVertices, mapping)){
				if (!used[e.toColour()][e.toNumber()] && uf.find(e.from()) != uf.find(e.to())){
					uf.union(e.from(), e.to());
					used[e.toColour()][e.toNumber()] = true;
					result.add(e);
					score += e.getScore();
				}
			}
    }


    return new AlgoResult(score, result);
  } */

	/*  public static void BnB(Vertex PMD, List<Edge> edgesList, int peakCount, int numberCount, int numberVertices){

    Edge currEdge, removedEdge;
    Edge[] edges = edgesList.toArray(new Edge[edgesList.size()]);
    Arrays.sort(edges);
    double sumOfScores = 0.0, currentScore = 0.0, possibleMaxScore = 0.0, currentMaxScore = 0.0;
    int[] colours = new int[peakCount];
    int[] next = new int[edges.length];
    boolean[][] used = new boolean[peakCount][numberCount];

    boolean fromOK, toOK, resultsToBeCleared = false;
    int scoredEdges = 0, nextInc = 0;
    LinkedList<Edge> result = new LinkedList<Edge>();

    for (int i = 0; i < colours.length; ++i){
      colours[i] = -1;
    }

    Edge[] currentEdges = new Edge[edges.length];
    int level = 0, currentEdgesPos = 0;

    // count maxScore
    for (int i = 0; i < edges.length; ++i){
      sumOfScores += edges[i].score;
    }

    while (level != -1 && currentMaxScore < sumOfScores){
      if (level != currentEdgesPos) System.out.println("Problem: level != currEdges");
      // calculate remainingScore
      possibleMaxScore = currentScore;
      scoredEdges = currentEdgesPos;
      // can't have more edges in result than are colours - 1
      for (int i = 0; i < edges.length && scoredEdges < peakCount-1; ++i){
        currEdge = edges[i];
        fromOK = colours[currEdge.fromColour()] == -1 || currEdge.fromNumber() == colours[currEdge.fromColour()];
        toOK = colours[currEdge.toColour()] == -1 || currEdge.toNumber() == colours[currEdge.toColour()];
        if (!used[currEdge.toColour()][currEdge.toNumber()] && fromOK && toOK) {
          possibleMaxScore += edges[i].score;
          ++scoredEdges;
        }
      }
      // if the remaining score can not reach the maxScore or
      // there are no more edges for this level, remove last edge from current edges
      if (next[level] >= edges.length || possibleMaxScore < currentMaxScore){
        //if (next[level] >= edges.length) System.out.println("cause array overrun");
        //if (possibleMaxScore < currentMaxScore) System.out.println("cause score trouble");
        --level;
        resultsToBeCleared = true;
        if (level != -1){
          --currentEdgesPos;
          removedEdge = currentEdges[currentEdgesPos];
          //System.out.println("Removed edge: "+removedEdge.toString());
          currentScore -= removedEdge.score;
          used[removedEdge.toColour()][removedEdge.toNumber()] = false;
          if (removedEdge.toHasSetColour) {
            colours[removedEdge.toColour()] = -1;
            removedEdge.toHasSetColour = false;
          }
          if (removedEdge.fromHasSetColour) {
            colours[removedEdge.fromColour()] = -1;
            removedEdge.fromHasSetColour = false;
          }
        }

      } else{
        // check whether edge is ok to add
        currEdge = edges[next[level]];
        fromOK = colours[currEdge.fromColour()] == -1 || currEdge.fromNumber() == colours[currEdge.fromColour()];
        toOK = colours[currEdge.toColour()] == -1 || currEdge.toNumber() == colours[currEdge.toColour()];
        if (!used[currEdge.toColour()][currEdge.toNumber()] && fromOK && toOK){
          currentScore += currEdge.score;
          if (colours[currEdge.fromColour()] == -1){
            currEdge.fromHasSetColour = true;
            colours[currEdge.fromColour()] = currEdge.fromNumber();
          }
          if (colours[currEdge.toColour()] == -1){
            currEdge.toHasSetColour = true;
            colours[currEdge.toColour()] = currEdge.toNumber();
          }
          used[currEdge.toColour()][currEdge.toNumber()] = true;
          currentEdges[currentEdgesPos] = currEdge;
          //System.out.println("Added edge: "+currEdge.toString());
          ++currentEdgesPos;
          ++next[level];
          ++level;
          if (level < next.length) next[level] = 0;
          if (currentScore > currentMaxScore){
            //System.out.println("Curr Max = "+currentScore);
            currentMaxScore = currentScore;
            if(resultsToBeCleared){ 
              result.clear();
              for (int i=0; i<currentEdgesPos; ++i){
                result.add(currentEdges[i]);
              }
              resultsToBeCleared = false;
            } else {
              result.add(currEdge);
            }
          }
        } else { // edge is not added, try next one
          ++next[level];
        }
      }
    }

  	PMD.setScore(PMD.getScore()+currentMaxScore);
  	PMD.setOptimalTree(result);
	}*/
}

class AlgoResult implements Comparable<AlgoResult> {
	public double score;
	public List<Edge> graph;
	public AlgoResult(double score, List<Edge> graph){
		this.score = score;
		this.graph = graph;
	}

	public int compareTo(AlgoResult another){
		if (score > another.score) return 1;
		if (score == another.score) return 0;
		return -1;    
	}
}

class TransitiveEdgeColor{
	public double score;
	public int from;
	public int to;

	public TransitiveEdgeColor(int from, int to, double score){
		this.from=from;
		this.to=to;	
		this.score=score;	
	}

	public int from(){
		return from;
	}

	public int to(){
		return to;
	}

	public String toString(){
		return from+" "+to+" "+score;
	}
}

class EdgeColor{
	public int from;
	public int to;

	public EdgeColor(int from, int to){
		this.from=from;
		this.to=to;	
	}

	public int from(){
		return from;
	}

	public int to(){
		return to;
	}

	public String toString(){
		return from+" "+to;
	}
}


