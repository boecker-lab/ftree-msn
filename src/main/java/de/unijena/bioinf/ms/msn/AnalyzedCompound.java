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
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Iterator;

/**
 * @author Kerstin Scheubert
 */
public class AnalyzedCompound{
  private String name;
  private Peak parentPeak;
  private List<Peak> peaks;  
  private List<Peak> unscoredPeaks=new ArrayList<Peak>();  

  public AnalyzedCompound(String name, Peak parentPeak, List<Peak> peaks){
    this.name = name;
    this.parentPeak = parentPeak;
    this.peaks = peaks;
  }

  public AnalyzedCompound(Compound target, Parameters param){
    name = target.getName();
    peaks = target.mergePeaks(param);    
    
    //adjust focused mass
    // DBE Check tests whether formula is already charged    
    if (target.getFormula() != null && Math.abs(target.getFormula().getMass() - target.getFocusedMass()) < 0.01 && target.getFormula().isNeutral()){
      target.setFocusedMass(target.getFocusedMass()+(target.getCharge()*Element.PROTON_MASS));
      //System.out.println("adjusted");
    }
    
    //remove isotopes from the merged spectrum
    // this was only for some weird data from Halle
    /*Spectrum mergedSpectrum = new Spectrum(0, peaks);
    if (param.isotopes || param.removeiso){
      mergedSpectrum.findIsotopePeaks(param);
    }*/
     
    parentPeak = null;    
    // look for a parent peak, closest to focused mass
    double diff = 0.0, minDiff = Double.POSITIVE_INFINITY, maxIntensity=0.0;    
    if (target.getFocusedMass() != 0.0){
      for (Peak p : peaks){
       /* diff = Math.abs(p.getMass()-target.getFocusedMass());
        if (diff < minDiff){
          minDiff = diff;
          parentPeak = p;
        }*/    	
    	diff = Math.abs(p.getMass()-target.getFocusedMass());    
        if (diff<=param.parentPeakError&&(parentPeak==null||parentPeak.getIntensity()<p.getIntensity())){
          minDiff = diff;
          parentPeak = p;          
        }    	
      }          
                  
      //TODO: change for absolute masses   
      double expectedParentMassNoH = 0.0;
	  if(param.input.getFormula()!=null)expectedParentMassNoH=target.getFormula().getMass();	  
	  double expectedParentMass = target.getCharge()*Element.PROTON_MASS+expectedParentMassNoH;      
	  if (parentPeak==null||(target.getFormula()!=null&&param.errorPpm<Math.abs(parentPeak.getMass()-expectedParentMassNoH)/parentPeak.getMass()*1e6&&param.errorPpm<Math.abs(parentPeak.getMass()-expectedParentMass)/parentPeak.getMass()*1e6)) {
        if (target.getMS1Peaks() != null && !target.getMS1Peaks().isEmpty()){
          parentPeak=null;
        } else if (target.getFormula() != null){
          parentPeak = new Peak(target.getFocusedMass(), 1, 0);
          parentPeak.setRelIntensity(parentPeak.getIntensity());
          peaks.add(parentPeak);
          Collections.sort(peaks);
          System.out.println("Parent peak generated");
        } else {
          System.out.println("No parent peak found");
          throw new RuntimeException("No parent peak found");
        }
      }
    } else { // use last peak if parent peak is unknown.
      parentPeak = peaks.get(peaks.size()-1);
    }

    List<Peak> ms1peaks = target.getMS1Peaks();
    if (ms1peaks != null && !ms1peaks.isEmpty()){
      if (parentPeak != null) {
        peaks.remove(parentPeak);
        // TODO: for better intensity normalisation, isotope anaylsis will crash!
        ms1peaks.get(0).setIntensity(parentPeak.getIntensity());
        ms1peaks.get(0).setHighestEnergy(parentPeak.getHighestEnergy());
        ms1peaks.get(0).setLowestEnergy(parentPeak.getLowestEnergy());
      } else {
        ms1peaks.get(0).setHighestEnergy(0);
        ms1peaks.get(0).setLowestEnergy(0);				
      }
      parentPeak = ms1peaks.get(0);
      peaks.add(parentPeak);
      for (Peak p : ms1peaks.subList(1, ms1peaks.size()))  {
        //p.setRelIntensity(p.getIntensity());
        parentPeak.addIsotope(p);
      }

      // erase all peaks > parentPeak
      for(Iterator<Peak> iter = peaks.iterator(); iter.hasNext();){
        if (iter.next().getMass() > parentPeak.getMass()){
          iter.remove();
        }
      }

    }

    if (param.parentiso) {
      Spectrum spectOfParent = target.getSpectra().get(parentPeak.getLowestEnergy());
      
      parentPeak.scoreIsotopeDist(param, spectOfParent, target.getCharge());
    }

    for (Iterator<Peak> iter = peaks.iterator(); iter.hasNext();){
      Peak p = iter.next();
      if (p != parentPeak && p.getIntensity()<param.intensityCutoff){
        iter.remove();
      }
    }
    
    //find all parentPeaks
    for(Spectrum s:target.getSpectra()){
    	if(s.getMSLevel()==2){
			s.setParentPeak(parentPeak);
		}
    	else{
	    	for(Peak p:peaks){	    		
    			//boolean isBestParentPeak=(Math.abs(s.getParentMass()-p.getMass())<param.parentPeakError)&&(s.getParentPeak()==null||Math.abs(s.getParentMass()-s.getParentPeak().getMass())>Math.abs(s.getParentMass()-p.getMass()));
    			boolean isBestParentPeak=(Math.abs(s.getParentMass()-p.getMass())<param.parentPeakError)&&(s.getParentPeak()==null||s.getParentPeak().getIntensity()<p.getIntensity());
				if(isBestParentPeak){
					s.setParentPeak(p);
				}	    
	    	}
    	}
    }     
  }

  public String getName(){
    return name;
  }

  public List<Molecule> getDecompositions(){
    return parentPeak.getDecompositions();
  }

  public List<Peak> getPeaks(){
    return peaks;
  }
  
  public List<Peak> getUnscoredPeaks(){
	    return unscoredPeaks;
}

  public void setName(String name){
    this.name = name;
  }

  public void setPeaks(List<Peak> peaks){
    this.peaks = peaks;
  }

  public Peak getParentPeak() {
    return parentPeak;
  }

  public void setParentPeak(Peak parentPeak) {
    this.parentPeak = parentPeak;
  }

  public Graph constructGraph(Molecule PMD, Parameters param, double parentmass) throws ElementNotFoundException{
    int peakCount = 0, decompCount = 0;
    boolean isParentPeak = false;
    List<Vertex> vertices = new ArrayList<Vertex>();
    List<Edge> edges = new ArrayList<Edge>();
    Vertex newVertex;

    //sort peaks in descending Intensity order 
    Collections.sort(peaks, Peak.relIntensityComparator());
    
    int peaksExact = getPeaks().size();
    if (PMD == null){
      peaksExact = Math.min(getPeaks().size(), param.peaksToScore);
    }

    ArrayList<Peak> scoredPeaks=new ArrayList<Peak>(); 
    scoredPeaks.add(parentPeak);
    if(!param.mostIntensivePeaks){
	    for(Spectrum s:param.input.getSpectra()){
	    	if(!scoredPeaks.contains(s.getParentPeak())){
	    		scoredPeaks.add(s.getParentPeak());    			
	    	}
	    }    
    }
       
    Iterator<Peak> pi=peaks.iterator();
    while(pi.hasNext()&&scoredPeaks.size()<peaksExact){
    	Peak nextPeak=pi.next();
    	if(!scoredPeaks.contains(nextPeak)){
    		scoredPeaks.add(nextPeak);    		
    	}
    }   
    
    Collections.sort(scoredPeaks,Peak.relIntensityComparator());
    scoredPeaks.remove(parentPeak);
    scoredPeaks.add(0,parentPeak);
    for(int i=scoredPeaks.size()-1;i>=Math.min(peaksExact,scoredPeaks.size());i--){
    	scoredPeaks.remove(i);
    }
    
    unscoredPeaks=new ArrayList<Peak>();    
    for(Peak p:peaks){
    	unscoredPeaks.add(p);
    }
    unscoredPeaks.removeAll(scoredPeaks);
    
    //System.out.println(getPeaks().subList(0, peaksExact))

    // create vertices  
//    Collections.sort(scoredPeaks,Peak.massComparator());
//    for(Peak p:scoredPeaks){
//    	System.out.println(p);
//    } 
    
    for (Peak peak : scoredPeaks){
      isParentPeak = (peak == parentPeak);
      decompCount = 0;
      for (Molecule decomposition : peak.getDecompositions()){
        if (PMD == null || decomposition.isSubsetOf(PMD)){
          boolean isSubsetOfParentPeakDecomposition=false;
          for(Molecule parentPeakDecomposition:parentPeak.getDecompositions()){
        	  if(decomposition.isSubsetOf(parentPeakDecomposition))isSubsetOfParentPeakDecomposition=true;        		 
          }
          if(isSubsetOfParentPeakDecomposition){
	          newVertex = new Vertex(decomposition, peak, peakCount, decompCount, isParentPeak);
	          peak.setColour(peakCount);
	          decomposition.setVertex(newVertex);
	          vertices.add(newVertex);
	          ++decompCount;
          }
        }
      }
      //if (decompCount != 0) {
        ++peakCount;
      //}
    }
  
    //adds a colour to all unscored peaks (only important for scoring msn spectra)
    Collections.sort(unscoredPeaks,Peak.relIntensityComparator());
    for(Peak p:unscoredPeaks){
    	p.setColour(peakCount);
    	++peakCount;
    }

    // create edges
    Edge newEdge;
    for (Vertex v : vertices){
      for (Vertex v2 : vertices){
        if (v2.getLabel().isTrueSubsetOf(v.getLabel())){
          //System.out.println("vertices "+v+" "+v2);
          newEdge = new Edge(v, v2, parentmass, param);
          //if (newEdge.score > 50 || v.isPMD()) {
          edges.add(newEdge);
          v.addOutEdge(newEdge);
          v2.addInEdge(newEdge);
          //}
        }
      }
    }

    /*to slow
    // erase all vertices not linked with any parent mass decomposition.
    if (PMD == null){
      boolean erasedSomething = true;
      Vertex v = null;
      while (erasedSomething) {
        erasedSomething = false;
        for (Iterator<Vertex> vIter = vertices.iterator(); vIter.hasNext(); ){
          v = vIter.next();
          if (!v.isPMD() && v.inEdges().size() == 0){
            v.delete();
            vIter.remove();
            erasedSomething = true;
          }
        }
      }
    }

    for (Iterator<Edge> eIter = edges.iterator(); eIter.hasNext();){
      if (eIter.next().isDeleted()){
        eIter.remove();
      }
    }
    
    */
    
    /*for (Vertex v : vertices){
    	v.rankOutEdges();
    }*/  
    if(!param.noTransitiveClosureScore||param.showTransitiveClosure){
      HashMap<Integer,HashMap<Integer,Double>> ms3score=new HashMap<Integer,HashMap<Integer,Double>>();
	  Collections.sort(peaks,Peak.massComparator());    
	  for(Peak cp:peaks){
		cp.setMS3Score(ms3score);		
		for(Spectrum s:cp.getSpectra()){			
			if(s.getMSLevel()>2&&s.getParentPeak()!=cp)addEntry(s.getParentPeak().getColour(), cp.getColour(), param.closureScore[0], ms3score);	
			for(Peak pp:peaks){
				if (pp.getMass()>cp.getMass()&&pp.getMass()<=s.getParentPeak().getMass()){
						if(!pp.getSpectra().contains(s)){
							addEntry(pp.getColour(), cp.getColour(), param.closureScore[1], ms3score);
						}
//						else{
//							for(Spectrum s2:pp.getSpectra()){
//								if(s2.getCollisionEnergy()>=s.getCollisionEnergy()&&!cp.getSpectra().contains(s2)){
//									addEntry(pp.getColour(), cp.getColour(), param.closureScore[2], ms3score);									
//								}
//							}
//						}
				}
				if (pp.getMass()<cp.getMass()&&!pp.getSpectra().contains(s)){
					boolean found=false;
					for(Spectrum s2:pp.getSpectra()){
						if(cp.getSpectra().contains(s2)&&s2.getCollisionEnergy()<=s.getCollisionEnergy()){
							found=true;									
						}
					}
					if(found)addEntry(cp.getColour(), pp.getColour(), param.closureScore[2], ms3score);
				}
			}
		}
	  }
    }
	return new Graph(vertices, edges);
  }

  public int countDecompositions() {
    int result = 0;
    for (Peak p : peaks){
      result += p.getDecompositions().size();
    }
    return result;
  }

  public void filterDecompToPMD(List<Molecule> PMDs) {
    for (Peak p : peaks){
      for (Iterator<Molecule> i = p.getDecompositions().iterator(); i.hasNext();){
        Molecule m = i.next();
        boolean delete = true;
        for (Molecule PMD : PMDs){
          if (m.isSubsetOf(PMD)){
            delete = false;
          }
        }
        if (delete){
          i.remove();
        }
      }
    }
  }
  
  private void addEntry(int parentcolour, int childrencolour, double score, HashMap<Integer,HashMap<Integer,Double>> ms3score){
	  if(ms3score.get(parentcolour)==null)ms3score.put(parentcolour,new HashMap<Integer,Double>());
	  if(ms3score.get(parentcolour).get(childrencolour)==null)ms3score.get(parentcolour).put(childrencolour,0.0);
	  ms3score.get(parentcolour).put(childrencolour,ms3score.get(parentcolour).get(childrencolour)+score);
  }

}
