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

import java.text.DecimalFormat;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.FileNotFoundException;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Locale;
import java.util.Set;
import java.util.regex.Pattern;

/**
 * @author Kerstin Scheubert
 */
public class MS2Analyzer{
	
	public static void main(String[] args) throws ElementNotFoundException {
//		String mol="NOSCAPINE";		
//		String inputString="D:\\MS3\\MM48_MS3\\"+mol+"\\"+mol+".msn -ar -g D:\\MS3\\MM48_MS3\\"+mol+"\\"+mol+".dot -mi -raw off";
//		args=inputString.split(" ");
		
		Locale.setDefault(Locale.US);
		long graphTime = 0, algoTime = 0;
		// Parse Arguments
		Parameters param = parseArgs(args);
		if (!param.success) return;
		
		//if (param.input.getFormula() == null) return;

		//    adds a perfect parent peak to the spectra   
		//    List<Peak> ms1peaks = new ArrayList<Peak>(1);
		//    ms1peaks.add(new Peak(param.input.getFormula().getMass()+param.input.getCharge()*Element.PROTON_MASS, 10000, 0));
		//    param.input.setMS1Peaks(ms1peaks);    

		// Start main processing
		Decomposer decomp = new Decomposer(param);
	    Edge.init(param, param.nlCombinations, param.et);

		System.out.println("Decomposer constructed");

		if (param.raw) param.input.calculateRaw();
		if (param.relative) param.input.calculateSNR();
		if (!param.raw && !param.relative) param.input.copyIntensities();
		param.input.normIntensities2(10000); // normalize Intensities in any case.		
		param.input.strictFilter(param.intensityCutoff);
		System.out.println();
		long startTime = System.nanoTime();
		if (param.isotopes || param.removeiso) {
			param.input.findIsotopePeaks(param);
		}
		if (param.parentiso){
			System.out.println("WARNING: MS1-Isotope calculation disturbed in AnalyzedCompound-Constructor");
		}
		param.input.decomposePeaks(decomp, param);
		param.input.filterDecompositions();
		param.input.scoreDecompositions(param);
		AnalyzedCompound anaCompound = new AnalyzedCompound(param.input, param);
		Peak parentPeak = anaCompound.getParentPeak();

		parentPeak.decompose(decomp, param, param.input.getCharge());
		//  the loop below removes all parentDecompositions except for the correct one.
		int c=0;
		for (Iterator<Molecule> iter = parentPeak.getDecompositions().iterator(); iter.hasNext();) {
			Molecule p = iter.next();
			// Use commmented if for degradation analysis (Applied Data) 
			//    if ((param.input.getFocusedMass() < param.input.getFormula().getMass() && p.isSubsetOf(param.input.getFormula())) ||
			//    (param.input.getFocusedMass() > param.input.getFormula().getMass() && param.input.getFormula().isSubsetOf(p))){
			if (!((param.input.getFormula()==null&&param.input.getFormula2()==null)||(param.input.getFormula()!=null&&p.isSubsetOf(param.input.getFormula()))||(param.input.getFormula2()!=null&&p.isSubsetOf(param.input.getFormula2())))){
				iter.remove();
			}			
			++c;
		}
		
		parentPeak.scoreDecompositions(param, param.input.getSpectra().get(0), param.input.getCharge());

		List<Molecule> PMDs = anaCompound.getDecompositions();
		Collections.sort(PMDs, new Comparator<Molecule>(){

			public int compare(Molecule m1, Molecule m2) {
				return (int) Math.signum(m2.getIsotopeScore() - m1.getIsotopeScore());
			}

		});
		//    int isopos=1;
		//    for (Molecule PMD : PMDs){
		//      Molecule correctHadj= new Molecule(param.input.getFormula());
		//      correctHadj.addElement("H", -1);
		//      if (PMD.equals(param.input.getFormula()) || PMD.equals(correctHadj)){
		//        System.out.println(anaCompound.getName()+" "+ PMD.getMass() +" isotopes correct at "+isopos+" of "+PMDs.size());
		//      }
		//      ++isopos;
		//    }
		if ((param.parentiso || param.isotopes) && !PMDs.isEmpty() && PMDs.get(0).getIsotopeScore() > -8){
			boolean zeros = true;
			for (int i = 0; i<5 && i<PMDs.size(); ++i){
				if (Math.abs(PMDs.get(i).getIsotopeScore())>0.1){
					zeros = false;
				}
			}
			if (!zeros){
				//        int PMDLimit = 5;
				//        if (PMDs.size() > PMDLimit) {
				//          System.out.println("Using isotopes as filter");
				//          List<Molecule> newPMDs = new ArrayList<Molecule>(PMDLimit);
				//          for (Molecule PMD : PMDs.subList(0, PMDLimit)) {
				//            newPMDs.add(PMD);
				//          }
				//          PMDs = newPMDs;
				//          System.out.println("new size of PMDs: " + PMDs.size());
				//          // remove the molecules which do not fit the remaining PMDs.
				//          anaCompound.filterDecompToPMD(PMDs);
				//        } 

				for (Molecule PMD : PMDs){
					PMD.addScore(5*PMD.getIsotopeScore()+param.isobonus);
				}
			} else { System.out.println("Isotopes disabled");}
		} else { System.out.println("Isotopes disabled");}


		// for counting remove all decompositions not subset of the corrcet formula
		//    for (Peak p : anaCompound.getPeaks()){
		//      for (Iterator<Molecule> iter = p.getDecompositions().iterator(); iter.hasNext();) {
		//        Molecule m = iter.next();
		//        if (!m.isSubsetOf(param.input.getFormula())){
		//          iter.remove();
		//        }
		//      }      
		//    }

		// write number ofdecompositions and exit    
		//    try {
		//      BufferedWriter decomps = new BufferedWriter(new FileWriter("NumberDecomps.txt", true));
		//      decomps.write("\""+anaCompound.getName()+"\", ");
		//      for (Peak p : anaCompound.getPeaks()){
		//                if (!p.equals(anaCompound.getParentPeak())) {
		//                  decomps.write(p.getDecompositions().size() + ", ");
		//                }                
		//      }
		//      decomps.newLine();
		//      //decomps.write(anaCompound.getName()+", "+anaCompound.getPeaks().size()+"\n");
		//      decomps.close();
		//      return;
		//    } catch (Exception e) {
		//      e.printStackTrace();
		//    }

		System.out.println("Parent Peak: "+parentPeak.getMass()+" Intensity: "+parentPeak.getIntensity()+" First appears in spectrum "+parentPeak.getLowestEnergy()+" Last appears in spectrum "+parentPeak.getHighestEnergy());    
		long decompTime = System.nanoTime() - startTime;
		double expectedParentMassNoH=0.0;
		if(param.input.getFormula()!=null)expectedParentMassNoH = param.input.getFormula().getMass();
		double expectedParentMass = param.input.getCharge()*Element.PROTON_MASS+expectedParentMassNoH;
		System.out.println("Exp. parent mass: "+expectedParentMass+" Accuracy: "+(Math.abs(parentPeak.getMass()-expectedParentMass)/parentPeak.getMass()*1e6)+" ppm");
		System.out.println("Exp. parent mass noH: "+expectedParentMassNoH+" Accuracy: "+(Math.abs(parentPeak.getMass()-expectedParentMassNoH)/parentPeak.getMass()*1e6)+" ppm");
		// Decomposition done

		int remainingDecompositions = anaCompound.countDecompositions(); //remaining Decomps after the filtering
		int colours = anaCompound.getPeaks().size();
		System.out.println("Number of parent mass decompositions: "+PMDs.size());
		System.out.println("Vertices: "+remainingDecompositions+ " Colours: "+colours);
		System.out.println("Starting main analysis (Note that it is a lot faster on the first colors)");

		List<Vertex> vertices = null;
		boolean MST = param.gheuristic || param.tdheuristic;
		if (!MST) {
			// Start graph construction, score decompositions
			startTime = System.nanoTime();
			// Construct Graph also scores decompositions/vertices
			Graph g = anaCompound.constructGraph(null, param, parentPeak.getMass());
			vertices = g.getVertices(); // null to construct the graph including all PMDs
			System.out.println("Edges: "+g.getEdges().size());
			graphTime = System.nanoTime() - startTime;
			//	    DecimalFormat format = new DecimalFormat("0.000"), shorty = new DecimalFormat("0.0");
			//      try {
			//        BufferedWriter fullgraph = new BufferedWriter(new FileWriter("fullgraph.dot"));
			//        fullgraph.write("digraph "+anaCompound.getName()+"{\n");
			//        // for better collision energy output
			//        int[] ces = new int[param.input.getSpectra().size()];
			//        int k = 0;
			//        for (Spectrum s : param.input.getSpectra()){
			//          ces[k] = s.getCollisionEnergy();
			//          ++k;
			//        }
			//        int i= 0;
			//        for (Peak p : anaCompound.getPeaks()){
			//          fullgraph.write("subgraph cluster"+i+" {\n");
			//          fullgraph.write("label = \""+format.format(p.getMass())+" Da Int.: "+shorty.format(p.getRelIntensity()/100)+"\\nCE: ");
			//          for (int j = p.getLowestEnergy(); j <= p.getHighestEnergy(); ++j){
			//            fullgraph.write(Integer.toString(ces[j]));
			//            if (j != p.getHighestEnergy()){
			//              fullgraph.write(", ");
			//            }
			//          }
			//          fullgraph.write(" eV\"\n");
			//          for (Vertex v : vertices){
			//            if (v.getColour() == i){
			//              fullgraph.write(v.dotString(ces, param.input.getCharge()));
			//              fullgraph.newLine();
			//            }
			//          }
			//          fullgraph.write("}\n");
			//          ++i;
			//        }
			//        for (Edge e : g.getEdges()){
			//          //if (e.getScore() > 10){
			//            fullgraph.write(e.dotString(ces, param.input.getCharge())+"\n");
			//          //}
			//        }
			//        fullgraph.write("}\n");
			//        fullgraph.close();
			//        return;
			//      } catch (IOException e) {
			//        e.printStackTrace();
			//      }

			// Graph constructed

			/*for (Peak p : anaCompound.getPeaks()){
			  System.out.println(p.getMass()+" "+p.getDecompositions());
      }*/
			for (Molecule PMD : PMDs){
				scorePMD(PMD, param);
			}
			startTime = System.nanoTime();
			colours = Math.min(colours, param.peaksToScore);
			colorcodingDP(vertices, colours);

			if (anaCompound.getPeaks().size() > param.peaksToScore) {				
				List<Peak> unscoredPeaks=anaCompound.getUnscoredPeaks();
				if(param.attachRemaining){
					for (Molecule PMD : PMDs) {
						Tree t = new Tree(PMD.getVertex().optimalTree(true));
						attachRemaining(t, unscoredPeaks, param);
					}
				}
			}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
			algoTime = System.nanoTime() - startTime;

		} else { // Split the graph, use MST algorithm

			if (param.gheuristic){
				System.out.println("Using greedy heuristic");
			} else if (param.tdheuristic) {
				System.out.println("Using top-down heuristic");
			} else {
				System.out.println("Using brute force algorithm");
			}
			int maxDecompCount = 0;
			for (Peak peak : anaCompound.getPeaks()){
				maxDecompCount= Math.max(maxDecompCount, peak.getDecompositions().size());
			}

			if (!param.gheuristic && !param.tdheuristic) System.out.print("PMDs done (of "+PMDs.size()+"): ");
			int i = 1;
			for (Molecule PMD : PMDs) {

				//System.out.print(i+". PMD: "+PMD.toString()+" ");
				if (!param.gheuristic && !param.tdheuristic) System.out.print(i+" ");
				// Start graph construction, score decompositions
				startTime = System.nanoTime();
				// Construct Graph also scores decompositions/vertices
				Graph g = anaCompound.constructGraph(PMD ,param, parentPeak.getMass());
				List<Edge> edges = g.getEdges();
				int[][] numberMapping = new int[colours][maxDecompCount];
				int[] amountOfColour = new int[colours];
				for (Peak peak : anaCompound.getPeaks()){
					for (Molecule decomposition : peak.getDecompositions()){
						if (decomposition.isSubsetOf(PMD)) {
							numberMapping[decomposition.getVertex().getColour()][amountOfColour[decomposition.getVertex().getColour()]] = decomposition.getVertex().getNumber();
							++amountOfColour[decomposition.getVertex().getColour()];
						}
					}
				}

				graphTime += System.nanoTime() - startTime;
				// Graph constructed
				scorePMD(PMD, param);

				startTime = System.nanoTime();
				if (param.gheuristic){
					Heuristics.greedy(PMD.getVertex(), edges, colours, maxDecompCount);
				} else if (param.tdheuristic) {
					Heuristics.topDown(PMD.getVertex(), colours, maxDecompCount);
				} else {
					Graph.scoreGraphMST(PMD.getVertex(), g, amountOfColour, numberMapping, maxDecompCount,param.noTransitiveClosureScore);
				}

				//Graph.BnB(PMD.getVertex(), edges, colours, maxDecomp, g.getVertices().size());
				algoTime += System.nanoTime() - startTime;
				++i;
			} // end for PMDs
			System.out.println();

		} // end else splitting

		Collections.sort(PMDs, Molecule.scoreComparator());	  	

		// Output
		DecimalFormat df = new DecimalFormat("0.000");
		// Informational output
		System.out.println();
		Molecule correct = param.input.getFormula();
		Molecule correct2 = null;
		if (correct != null) {
			correct2 = new Molecule(correct);
			correct2.removeElement("H");
		}

		Molecule correct3 = param.input.getFormula2();
		Molecule correct4 = null;
		if (correct3 != null) {
			correct4 = new Molecule(correct);
			correct4.removeElement("H");
		}

		int correctPos = 0;
		if (!PMDs.isEmpty()) System.out.println("Best Decomposition: "+PMDs.get(0).toString()+" Score: "+df.format(PMDs.get(0).getScore()));
		if (param.input.getFormula() != null){
			int i = 1;
			for (Molecule PMD : PMDs){
				if (PMD.equals(correct) || PMD.equals(correct2) || PMD.equals(correct3) || PMD.equals(correct4)){
					correctPos = i;
					System.out.println(anaCompound.getName()+" "+ PMD.getMass() +" correct at "+i+" of "+PMDs.size());
					System.out.println("Correct decomposition: "+PMD.toString()+" Score: "+df.format(PMD.getScore()));
				}
				++i;        
			}
		}
		System.out.println("Decomposition time: "+decompTime/1e6+"ms Preprocessing time: "+graphTime/1e6+"ms\nAlgorithm time: "+algoTime/1e6+"ms");

		try {
			BufferedWriter scoref = new BufferedWriter(new FileWriter("scores.csv", true));
			BufferedWriter isoscoref = new BufferedWriter(new FileWriter("isoscores.csv", true));
			scoref.write(anaCompound.getName());
			isoscoref.write("\""+anaCompound.getName()+"\"");
			for (int i = 0; i<5 && i<PMDs.size(); ++i){
				Molecule m = PMDs.get(i);
				scoref.write(", "+m.getScore());
				isoscoref.write(", "+m.getIsotopeScore());
			}
			scoref.newLine();
			isoscoref.newLine();
			scoref.close();
			isoscoref.close();
		} catch (IOException e) {
			System.out.println("Could not write to table files");
		}

		// for better collision energy output
		int[] ces = new int[param.input.getSpectra().size()];
		int k = 0;
		for (Spectrum s : param.input.getSpectra()){
			ces[k] = s.getCollisionEnergy();
			++k;
		}

		// Output in file
		try {
			BufferedWriter outfile = new BufferedWriter(new FileWriter(param.outfile));
			BufferedWriter graphfile = null;
			outfile.write("Analysed spectra of "+param.input.getName()+"\n");
			outfile.newLine();
			outfile.write("Parent Peak: "+parentPeak.getMass()+" Da Intensity: "+parentPeak.getIntensity()+" First appears in spectrum "+parentPeak.getLowestEnergy()+" Last appears in spectrum "+parentPeak.getHighestEnergy());    
			outfile.newLine();
			outfile.write("Number of parent mass decompositions: "+PMDs.size());
			outfile.newLine();
			outfile.write("Vertices: "+remainingDecompositions+ " Colours: "+colours);
			outfile.newLine();
			outfile.newLine();
			param.write(outfile);
			outfile.newLine();
			outfile.newLine();

			if (!PMDs.isEmpty()){
				outfile.write("Best Decomposition: "+PMDs.get(0).toString()+" Score: "+df.format(PMDs.get(0).getScore()));
				outfile.newLine();
			}
			if (correctPos > 0){
				outfile.write("Correct at position: "+correctPos+" of "+PMDs.size());
				outfile.newLine();
				outfile.write("Correct decomposition: "+PMDs.get(correctPos-1).toString()+" Score: "+df.format(PMDs.get(correctPos-1).getScore()));
				outfile.newLine();
			}
			outfile.write("Decomposition time: "+decompTime/1e6+"ms Preprocessing time: "+graphTime/1e6+"ms\nAlgorithm time: "+algoTime/1e6+"ms");
			outfile.newLine();
			outfile.newLine();

			outfile.write(param.detailedOutput+" best sum formulas with their fragmentation trees:\n");
			for (int i = 0; i < param.detailedOutput && i < PMDs.size(); ++i){
				Molecule PMD = PMDs.get(i);
				outfile.write((i+1)+") "+ PMD.toString()+" "+df.format(PMD.getScore())+"\n");
				if (param.writeGraph) {
					String corrString = (PMD.equals(correct) || PMD.equals(correct2))?" correct":"";
					graphfile = new BufferedWriter(new FileWriter(param.graphfile+(i+1)+".dot"));
					graphfile.write("digraph "+PMD.toString()+" {\n");
					graphfile.write("label=\""+PMD.toString()+" "+df.format(PMD.getScore())+" "+(i+1)+"th suggest for "+anaCompound.getName()+corrString+"\"\n");
					graphfile.write("ranksep=3\n");
				}
				List<Edge> tree = PMD.getVertex().optimalTree(true);
				if (tree == null) System.out.println("tree = null");

				if (param.writeGraph) {
					Set<Vertex> treeVertices = new HashSet<Vertex>(tree.size() + 1);
					for (Edge e : tree) {
						treeVertices.add(e.from());
						treeVertices.add(e.to());
					}
					for (Vertex v : treeVertices) {
						graphfile.write(v.dotString(ces, param.input.getCharge()) + "\n");
					}
				}				

				for (Edge e : tree){
					outfile.write("\t"+e.toString()+"\n");
					if (param.writeGraph) graphfile.write(e.dotString(ces, param.input.getCharge())+"\n");
				}
				if (param.writeGraph){
					graphfile.write("}\n");
					graphfile.write("// Max Edge Score: "+Edge.getMaxEdgeScore(param));
					graphfile.close();
				}
				outfile.newLine();
			}
			if (correctPos > param.detailedOutput){
				Molecule PMD = PMDs.get(correctPos-1);
				outfile.write(correctPos+") "+ PMD.toString()+" "+df.format(PMD.getScore())+"\n");
				if (param.writeGraph) {
					String corrString = (PMD.equals(correct) || PMD.equals(correct2))?" correct":"";
					graphfile = new BufferedWriter(new FileWriter(param.graphfile+correctPos+".dot"));
					graphfile.write("digraph "+PMD.toString()+" {\n");
					graphfile.write("label=\""+PMD.toString()+" "+df.format(PMD.getScore())+" "+correctPos+"th suggest for "+anaCompound.getName()+corrString+"\"\n");
					graphfile.write("ranksep=3\n");
				}
				List<Edge> tree = PMD.getVertex().optimalTree(true);
				if (tree == null) System.out.println("tree = null");
				for (Edge e : tree){
					outfile.write("\t"+e.toString()+"\n");
					if (param.writeGraph) graphfile.write(e.dotString(ces, param.input.getCharge())+"\n");
				}
				if (param.writeGraph){
					graphfile.write("}\n");
					graphfile.close();
				}
				outfile.newLine();        
			}
			outfile.newLine();
			outfile.write("Complete list of parentmass decompositions with scores: \n");
			int i = 1;
			for (Molecule PMD : PMDs){
				outfile.write(i+") "+PMD.toString()+" "+df.format(PMD.getScore())+"\n");
				++i;
			}
			outfile.close();	    
			System.out.println("Detailed output written to "+param.outfile);

			/*BufferedWriter mergedspecFile = new BufferedWriter(new FileWriter(anaCompound.getName()+"Merged.ms"));
      DecimalFormat none = new DecimalFormat("0");
      for (Peak p : anaCompound.getPeaks()){
        mergedspecFile.write(df.format(p.getMass()));
        mergedspecFile.write(" ");
        mergedspecFile.write(none.format(p.getRelIntensity()));
        mergedspecFile.newLine();
      }
      mergedspecFile.close();*/
		}
		catch (IOException e){
			System.err.println("Could not write to output file: "+e.toString());
		} // */
	} // end function

	private static void colorcodingDP(List<Vertex> vertices, int numberOfColours){

		// Inverse Topologically sort vertices
		Collections.sort(vertices);

		for (Vertex v : vertices){
			v.initSets(numberOfColours);
		}

		System.out.print("Colours processed: 1 ");
		for (int sizeOfSet = 2; sizeOfSet <= numberOfColours; ++sizeOfSet){
			System.out.print(sizeOfSet+" ");
			for (Vertex v : vertices){
				// First part of update: use new edges
				for (Edge oe : v.outEdges()){
					Vertex neighbour = oe.to();
					// walk over all sets with size (sizeOfSet-1)
					DoubleIterator vIter = neighbour.getSets().values(sizeOfSet-1);
					for (LongIterator kIter = neighbour.getSets().keys(sizeOfSet-1); kIter.hasNext();){
						long newSet = BitSet.add(v.getColour(), kIter.next());
						v.getSets().updateSet(newSet, sizeOfSet, vIter.next() + oe.getScore());
					}
				}

				// Second part of update
				// cardinality(S1) + cardinality(S2) = sizeOfSet+1, cause colout(v) is twice in the sets.
				for (int sizeOfS1 = 2; sizeOfS1 <= (sizeOfSet+1)/2.0; ++sizeOfS1){
					DoubleIterator vIter = v.getSets().values(sizeOfS1);
					for (LongIterator kIter = v.getSets().keys(sizeOfS1); kIter.hasNext();){
						long S1 = kIter.next();
						double scoreS1 = vIter.next();
						DoubleIterator vIter2 = v.getSets().values(sizeOfSet+1-sizeOfS1);
						for (LongIterator kIter2 = v.getSets().keys(sizeOfSet+1-sizeOfS1); kIter2.hasNext();){
							long S2 = kIter2.next();
							double scoreS2 = vIter2.next();
							if (BitSet.size(BitSet.intersect(S1, S2)) == 1){
								v.getSets().updateSet(BitSet.union(S1, S2), sizeOfSet, scoreS1 + scoreS2);
							}
						}
					}
				}
			}
		}
		System.out.println("");

		// Assign scores to the labels of the sources, namely the PMDs
		for (Vertex source : vertices) {
			if (source.inEdges().isEmpty()){
				source.addScore(source.getSets().maxScore());
			}
		}

	}

	static void attachRemaining(Tree tree, List<Peak> peaks, Parameters param){
//		for (Edge e : tree.getEdges()){
//			System.out.println(e.dotString());
//		}

		if (tree.getEdges().isEmpty() || peaks.isEmpty()){
			return;
		}
		//sort peaks in descending Intensity order 
		Collections.sort(peaks, Peak.relIntensityComparator());
		double parentmass = tree.getOrgRoot().getLabel().getMass(), sumGains = 0.0;
		for (Peak p : peaks){
			double maxScoreGain = 0.0, scoreGain = 0.0, diff=0.0;
			Vertex maxDecompV = null, maxParent = null, vertexM = null;
			Edge newEdge = null, newEdgeToChild = null;
			for (Molecule m : p.getDecompositions()){
				vertexM = new Vertex(m, p, -1, -1, false);
				for (Vertex v : tree.getVertices()){
					if (m.isTrueSubsetOf(v.getLabel())) {
						newEdge = new Edge(v, vertexM, parentmass, param);
						scoreGain = newEdge.getScore();
						for (Edge oe : v.outEdges()) {
							if (oe.to().getLabel().isTrueSubsetOf(m)) {
								newEdgeToChild = new Edge(vertexM, oe.to(), parentmass, param);
								diff = newEdgeToChild.getScore() - oe.getScore();
								if (diff > 0) {
									scoreGain += diff;
								}
							}              
						}
						if (scoreGain > maxScoreGain) {
							maxScoreGain = scoreGain;
							maxDecompV = vertexM;
							maxParent = v;
						}
					}          
				}
			}

			if (maxScoreGain > 0){
				newEdge = new Edge(maxParent, maxDecompV, parentmass, param);
				try {
					tree.addEdge(newEdge);
					sumGains += maxScoreGain;//newEdge.getScore();
					for(Edge oe: maxParent.outEdges()){
						if (oe.to().getLabel().isTrueSubsetOf(maxDecompV.getLabel())) {
							newEdgeToChild = new Edge(maxDecompV, oe.to(), parentmass, param);
							diff = newEdgeToChild.getScore() - oe.getScore();
							if (diff > 0) {
								tree.addEdgeRemoveCycle(newEdgeToChild);
								//                
								//                // This probably never happens: Remove node if this edge was the last outgoing edge, and node gives a negative score
								//                System.out.println("Vacant edge "+oe.dotString());
								//                System.out.println(oe.to().outEdges());
								//                if (oe.to().outEdges().isEmpty() && !oe.from().inEdges().isEmpty()){
								//                  Edge removalCand = oe.from().inEdges().get(0);
								//                  System.out.println("removable edge "+removalCand.dotString());
								//                  if (removalCand.getScore() <= 0){
								//                    removalCand.to().inEdges().remove(removalCand);
								//                    removalCand.from().outEdges().remove(removalCand);
								//                    tree.getEdges().remove(removalCand);
								//                  }
								//                }
								//sumGains += diff;
							}
						}          
					}
				} catch (RuntimeException e) {
					// noop
				}        
			}
		}

		// Do/Do not include the attached edges to overall score!
		//tree.getOrgRoot().addScore(sumGains);
		tree.getOrgRoot().setOptimalTree(tree.getEdges());

	}

	static Parameters parseArgs(String[] args){
		Parameters res = new Parameters();
		int i = 0;
		boolean outfileSet = false;
		BufferedReader sFile = null;

		while (i < args.length){
			if (args[i].equals("-h") || args[i].equals("--help")) {
				System.out.println("\nSyntax: MS2Analyzer [options] inputfile\n");
				System.out.println("Available Options:");
				System.out.println("Start every new option with a dash (-), combined options might be interpreted as another option.\n");
				System.out.println("-a <filename> -et <filename>");
				System.out.println("\tRead alphabet masses and valences from specified file.\n\tIf none is given CHNOPS is used.\n");

				System.out.println("-e <double>");
				System.out.println("\tSpecify the mass deviation in ppm. Default: 20 ppm\n");

				System.out.println("-ea <double>");
				System.out.println("\tSpecify the absolute mass deviation in mDa. Default: 1 mDa\n");

				System.out.println("-p <double>");
				System.out.println("\tSpecify the precision of the decomposition. Default: 1e-5\n");

				System.out.println("-raw on|off");
				System.out.println("\tEnable/disable calculation of raw intensities based on the TIC\n\tof the spectra. Enabled by default\n");

				System.out.println("-rel on|off");
				System.out.println("\tEnable/Disable smoothing of peak intensities by scoring the most intense\n\tpeak 1, the second intense peak");
				System.out.println("\t(#peaks in spectrum-1)/#peaks in spectrum and so on.\n");

				System.out.println("-m <double>");
				System.out.println("\tSpecify how close peaks have to be, to be merged together to one peak.\n");

				System.out.println("-mi");
				System.out.println("\tIf two peaks are merged, keep the mass of the most intense peak.");
				System.out.println("\tBy default the average mass is taken.");

				//System.out.println("-s");
				//System.out.println("\tDo not use the peak with highest mass as parent peak, but scan\n\tfor more intense peaks in the vicinity.\n");

				System.out.println("-md <double>");
				System.out.println("\tSpecify the multiples of standard deviation, which are\n\tregarded to be inside the mass deviation error.");
				System.out.println("\tThis is used to score the fragment decompositions.");
				System.out.println("\tDefault: 3, that is deviation is mass error divided by 3\n");

				System.out.println("-n <filename>");
				System.out.println("\tSpecify a file containing sum formulas of neutral losses,\n\twhich are likely to be lost by your compound.");
				System.out.println("\tA decent list for natural compounds is used by default.\n");

				System.out.println("-nc <int>");
				System.out.println("\tCombinations of how many likely neutral losses are\n\ttreated as likely neutral losses, too.\n\tDefault: 3");

				System.out.println("-hc [average deviation]");
				System.out.println("\tScore using Hydrogen to Carbon ratio.");
				System.out.println("\tThe first argument is the average of the underlying normal distribution, the second the standard deviation.");
				System.out.println("\tDefaults (used if values are omited): Average = 1.44 Deviation = 0.50\n");

				System.out.println("-he [average deviation]");
				System.out.println("\tScore using Hetero atom to Carbon ratio.");
				System.out.println("\tThe first argument is the average of the normal distribution, the second the standard deviation.");
				System.out.println("\tDefaults (used if values are omited): Average = 0.58 Deviation = 0.56\n");

				System.out.println("-d [average deviation]");
				System.out.println("\tScore using DBE distribution.");
				System.out.println("\tThe first argument is the average of the normal distribution, the second the standard deviation.");
				System.out.println("\tDefaults (used if values are omited): Average = 6.15 Deviation = 4.54\n");

				System.out.println("-gh");
				System.out.println("\tUse the greedy heuristic to calculate the fragmentation trees.");

				System.out.println("-th");
				System.out.println("\tUse the top-down heuristic to calculate the fragmentation trees.");

				System.out.println("-i");
				System.out.println("\tCheck spectrum for isotope peaks and use isotope pattern scoring.");        
				System.out.println("-pi");
				System.out.println("\tUse isotope pattern scoring of parent peak.");        
				System.out.println("-ri");
				System.out.println("\tRemove isotopic peaks from the spectrum.");        

				System.out.println("-ie <double>");
				System.out.println("\tSpecify intensity error for isotope scoring in percent.");				

				System.out.println("-t <double>");
				System.out.println("\tIntensity threshold/cutoff");								

				System.out.println("-f <int>");
				System.out.println("\tNumber of fragment trees to be printed in the output file. Default: 10\n");

				System.out.println("-g <filename>");
				System.out.println("\tWrite the fragment trees as dot-files. The argument is used as filename-prefix.\n");

				System.out.println("-o <filename>");
				System.out.println("\tName of the output file. Default: <inputfile>.out\n");


				return res;
			}else if (args[i].equals("-a") || args[i].equals("-et")){
				++i;
				if (i == args.length || args[i].charAt(0) == '-'){
					System.err.println("Option \""+args[i-1]+"\" requires an argument.");
					return res;
				}
				try {
					BufferedReader elFile = new BufferedReader(new FileReader(args[i]));
					res.et = readElementTable(elFile);
					elFile.close();
				} catch (FileNotFoundException fnfe){
					System.err.println("Element file \""+args[i]+"\" not found");
					return res;
				} catch (IOException e) {
					System.err.println("Error reading element file "+e.getMessage());
					return res;
				}
			}
			else if (args[i].equals("-e")){
				++i;
				if (i == args.length || args[i].charAt(0) == '-'){
					System.err.println("Option \""+args[i-1]+"\" requires an argument.");
					return res;
				}
				try { res.errorPpm = Double.parseDouble(args[i]); }
				catch (NumberFormatException e) {
					System.err.println("Could not parse \""+args[i]+"\", is it a numeric value?");
					return res;
				}
				if (res.errorPpm < 0) {
					System.err.println("Mass error needs to be positive");
					return res;
				}
			}
			else if (args[i].equals("-ae")){
				++i;
				if (i == args.length || args[i].charAt(0) == '-'){
					System.err.println("Option \""+args[i-1]+"\" requires an argument.");
					return res;
				}
				try { res.absError = Double.parseDouble(args[i]); }
				catch (NumberFormatException e) {
					System.err.println("Could not parse \""+args[i]+"\", is it a numeric value?");
					return res;
				}
				if (res.absError < 0) {
					System.err.println("Mass error needs to be positive");
					return res;
				}
			}
			else if (args[i].equals("-p")){
				++i;
				if (i == args.length || args[i].charAt(0) == '-'){
					System.err.println("Option \""+args[i-1]+"\" requires an argument.");
					return res;
				}
				try { res.precision = Double.parseDouble(args[i]); }
				catch (NumberFormatException e) {
					System.err.println("Could not parse \""+args[i]+"\", is it a numeric value?");
					return res;
				}
				if (res.precision < 0 || res.precision > 1) {
					System.err.println("Precision needs to be between 0 and 1");
					return res;
				}
			}
			else if (args[i].equals("-raw")){
				++i;
				if (i < args.length && args[i].equals("on")){
					res.raw = true;
				} else if (i < args.length && args[i].equals("off")) {
					res.raw = false;
				} else {
					System.err.println("Please specify \"on\" or \"off\" after -raw");
					return res;
				}
			}
			else if (args[i].equals("-rel")){
				++i;
				if (i < args.length && args[i].equals("on")){
					res.relative = true;
				} else if (i < args.length && args[i].equals("off")) {
					res.relative = false;
				} else {
					System.err.println("Please specify \"on\" or \"off\" after -rel");
					return res;
				}
			}
			else if (args[i].equals("-gh")){
				res.gheuristic = true;
			}
			else if (args[i].equals("-th")){
				res.tdheuristic = true;
			}
			else if (args[i].equals("-i")){
				res.isotopes = true;
			}
			else if (args[i].equals("-pi")){
				res.parentiso = true;
			}
			else if (args[i].equals("-ri")){
				res.removeiso = true;
			}
			else if (args[i].equals("-ie")){
				++i;
				if (i == args.length || args[i].charAt(0) == '-'){
					System.err.println("Option \""+args[i-1]+"\" requires an argument.");
					return res;
				}
				try { res.intError = Double.parseDouble(args[i]); }
				catch (NumberFormatException e) {
					System.err.println("Could not parse \""+args[i]+"\", is it a numeric value?");
					return res;
				}
				if (res.intError < 0 || res.intError > 100) {
					System.err.println("Intensity error should be between 0 and 100.");
					return res;
				}
			}
			else if (args[i].equals("-m")){
				++i;
				if (i == args.length || args[i].charAt(0) == '-'){
					System.err.println("Option \""+args[i-1]+"\" requires an argument.");
					return res;
				}
				try { res.mergeThreshold = Double.parseDouble(args[i]); }
				catch (NumberFormatException e) {
					System.err.println("Could not parse \""+args[i]+"\", is it a numeric value?");
					return res;
				}
				if (res.mergeThreshold < 0) {
					System.err.println("Threshold needs to be positive");
					return res;
				}
			}
			else if (args[i].equals("-mi")){
				res.intenseMerge = true;
			}
			else if (args[i].equals("-ar")){
				res.attachRemaining = true;
			}
			/*else if (args[i].equals("-s")){
				res.scan = true;
			}*/
			else if (args[i].equals("-md")){
				++i;
				if (i == args.length || args[i].charAt(0) == '-'){
					System.err.println("Option \""+args[i-1]+"\" requires an argument.");
					return res;
				}
				try { res.massDeviationPenalty = Double.parseDouble(args[i]); }
				catch (NumberFormatException e) {
					System.err.println("Could not parse \""+args[i]+"\", is it a numeric value?");
					return res;
				}
				if (res.massDeviationPenalty < 0) {
					System.err.println("Mass deviation penalty needs to be positive");
					return res;
				}			
			}
			else if (args[i].equals("-n")){
				++i;
				if (i == args.length || args[i].charAt(0) == '-'){
					System.err.println("Option \""+args[i-1]+"\" requires an argument.");
					return res;
				}
				try {
					BufferedReader nFile = new BufferedReader(new FileReader(args[i]));
					res.neutralLossList = readList(nFile);
					nFile.close();
				} catch (FileNotFoundException fnfe){
					System.err.println("Neutral loss file \""+args[i]+"\" not found");
					return res;
				} catch (IOException e) {
					System.err.println("Error reading neutral loss file "+e.getMessage());
					return res;
				}
			}
			else if (args[i].equals("-hc")){ // Complicated Structure: Read up to two numerical values, but if none follows, do not abort
				res.hcRatioScoring = true;
				if ((i+1) != args.length){
					try {
						res.hcAverage = Double.parseDouble(args[i+1]);
						++i;
						if ((i+1) != args.length){
							res.hcDev = Double.parseDouble(args[i+1]);
							++i;
						}								
					}	catch (NumberFormatException e) {}
				}
			}
			else if (args[i].equals("-he")){ // Complicated Structure: Read up to two numerical values, but if none follows, do not abort
				res.heteroScoring = true;
				if ((i+1) != args.length){
					try {
						res.heteroAverage = Double.parseDouble(args[i+1]);
						++i;
						if ((i+1) != args.length){
							res.heteroDev = Double.parseDouble(args[i+1]);
							++i;
						}								
					}	catch (NumberFormatException e) {}
				}
			}
			else if (args[i].equals("-d")){ // Complicated Structure: Read up to two numerical values, but if none follows, do not abort
				res.DBEDependentScoring = true;
				if ((i+1) != args.length){
					try {
						res.DBEAverage = Double.parseDouble(args[i+1]);
						++i;
						if ((i+1) != args.length){
							res.DBEStdDev = Double.parseDouble(args[i+1]);
							++i;
						}								
					}	catch (NumberFormatException e) {}
				}
			}
			else if (args[i].equals("-f")){
				++i;
				if (i == args.length || args[i].charAt(0) == '-'){
					System.err.println("Option \""+args[i-1]+"\" requires an argument.");
					return res;
				}
				try { res.detailedOutput = Integer.parseInt(args[i]); }
				catch (NumberFormatException e) {
					System.err.println("Could not parse \""+args[i]+"\", is it an integer?");
					return res;
				}
			}
			else if (args[i].equals("-nc")){
				++i;
				if (i == args.length || args[i].charAt(0) == '-'){
					System.err.println("Option \""+args[i-1]+"\" requires an argument.");
					return res;
				}
				try { res.nlCombinations = Integer.parseInt(args[i]); }
				catch (NumberFormatException e) {
					System.err.println("Could not parse \""+args[i]+"\", is it an integer?");
					return res;
				}
			}
			else if (args[i].equals("-t")){
				++i;
				if (i == args.length || args[i].charAt(0) == '-'){
					System.err.println("Option \""+args[i-1]+"\" requires an argument.");
					return res;
				}
				try { res.intensityCutoff = Double.parseDouble(args[i]); }
				catch (NumberFormatException e) {
					System.err.println("Could not parse \""+args[i]+"\", is it a number?");
					return res;
				}
			}
			else if (args[i].equals("-g")){
				++i;
				if (i == args.length || args[i].charAt(0) == '-'){
					System.err.println("Option \""+args[i-1]+"\" requires an argument.");
					return res;
				}
				res.writeGraph = true;
				res.graphfile = args[i];
				if (res.graphfile.endsWith(".dot")) res.graphfile = res.graphfile.substring(0,res.graphfile.length()-4);
			}
			else if (args[i].equals("-o")){
				++i;
				if (i == args.length || args[i].charAt(0) == '-'){
					System.err.println("Option \""+args[i-1]+"\" requires an argument.");
					return res;
				}
				res.outfile = args[i];
				outfileSet = true;
			} else { // take it as input file
				try {
					sFile = new BufferedReader(new FileReader(args[i]));
				} catch (FileNotFoundException fnfe){
					System.err.println("Input file \""+args[i]+"\" not found");
					return res;
				}
				if (!outfileSet) res.outfile = args[i]+".out";
			}
			++i;
		}

		if (res.tdheuristic && res.gheuristic) {
			System.err.println("May only use one heuristic. Please remove either \"-th\" or \"-gh\" option.");
			return res;
		}

		if (sFile != null){
			try {
				res.input = readSpectra(sFile, res);
				sFile.close();				
				res.success = true;
			} catch (IOException e) {
				System.err.println("Error reading input file "+e.toString());
				return res;
			}
		} else {
			System.err.println("No input file given, call 'MS2Analyzer -h' for a list of options");		
		}
		return res;
	}

	static ElementTable readElementTable(BufferedReader file) throws IOException{
		ElementTable res = new ElementTable();
		double weight;
		int valence, i = 1;
		for (String line = file.readLine(); line!=null; line = file.readLine()) {
			line = line.trim();
			if (!line.equals("")) {
				String[] values = line.split("\\s+");
				if (values.length < 3) throw new IOException("in line "+i+". Format is \"elementname weight valence\"");
				try { weight = Double.parseDouble(values[1]); }
				catch (NumberFormatException e) { throw new IOException("in line "+i+". String \""+values[1]+"\" could not be parsed to double"); }
				try { valence = Integer.parseInt(values[2]); }
				catch (NumberFormatException e) { throw new IOException("in line "+i+". String \""+values[2]+"\" could not be parsed to integer"); }
				res.add(new Element(values[0], weight, valence));
			}
			++i;
		}
		return res;
	}

	static String readList(BufferedReader file) throws IOException{
		String res = "";
		for (String line = file.readLine(); line!=null; line = file.readLine()){
			res += line+" ";
		}
		return res;
	}

	static Compound readSpectra(BufferedReader file, Parameters param) throws IOException{
		String name = "", str;
		Molecule formula = null, formula2 = null;
		int mode = 0, collisionEnergy = 0;
		double tic = 0.0, focusedMass = 0.0, ms1ParentMass = 0.0, nextRT = 0.0, currRT=0.0;
		List<Spectrum> spectra = new ArrayList<Spectrum>();
		List<Peak> peaks = null, ms1peaks = new ArrayList<Peak>(5);
		Pattern peakPattern = Pattern.compile("\\d+\\.?\\d*\\s+\\d+\\.?\\d*E?\\d*");
		int i = 1;
		for (String line = file.readLine(); line!=null; line = file.readLine()){
			line = line.trim();
			if (line.contains("compound")) {
				name = line.substring(line.indexOf("compound")+9);
			}
			else if (line.contains("formula")) {
				str = line.substring(line.indexOf("formula")+8);
				try { 
					if (formula == null) {
						formula = new Molecule(str, param.et);
					} else {
						formula2 = new Molecule(str, param.et);
					}
				}
				catch (ElementNotFoundException e){
					throw new IOException("in line "+i+". "+e.getMessage());				
				}
				catch (Exception e) { 
					throw new IOException("in line "+i+". \""+str+"\" is not a valid sum formula");
				}
			}
			else if (line.contains("charge")){ 
				str = line.substring(line.indexOf("charge")+7);
				try { mode = Integer.parseInt(str); }
				catch (NumberFormatException e) { throw new IOException("in line "+i+". String \""+str+"\" could not be parsed to integer.");}
			}
			else if (line.contains("parentmass")){ 
				if (line.contains("Da")){
					str = line.substring(line.indexOf("parentmass")+11, line.indexOf("Da")-1);
				} else {
					str = line.substring(line.indexOf("parentmass")+11);          
				}
				try { focusedMass = Math.max(Double.parseDouble(str),focusedMass); }
				catch (NumberFormatException e) { throw new IOException("in line "+i+". String \""+str+"\" could not be parsed to double.");}
			}			
			else if (line.contains("collision")){
				if (peaks != null && peaks.size() > 0){
					if (tic == 0.0 && param.raw) throw new IOException("in line "+i+". TIC of previous spectrum undefined, cannot calculate raw intensities.\nUse \"-nraw\" switch or correct input file.");
					spectra.add(new Spectrum(collisionEnergy, tic, peaks, currRT));
				}
				currRT = nextRT;
				int begin = 11, end = line.length();
				if (line.contains("CE")) { begin = line.indexOf("CE")+3; }
				if (line.contains("eV")) { end = line.indexOf("eV")-1; }
				else if (line.contains("//")) { end = line.indexOf("//")-1; }
				str = line.substring(begin, end);
				str.trim();
				try { collisionEnergy = Math.round(Float.parseFloat(str)); }
				catch (NumberFormatException e) { throw new IOException("in line "+i+". String \""+str+"\" could not be parsed to integer.");}
				tic = 0.0;
				peaks = new ArrayList<Peak>();
			}
			else if (line.contains("tic")){
				str = line.substring(line.indexOf("tic")+4);
				try { tic = Double.parseDouble(str); }
				catch (NumberFormatException e) { throw new IOException("in line "+i+". String \""+str+"\" could not be parsed to double.");}       
			}
			else if (line.contains("retention")){
				str = line.substring(line.indexOf("retention")+10);
				try { nextRT = Double.parseDouble(str); }
				catch (NumberFormatException e) { throw new IOException("in line "+i+". String \""+str+"\" could not be parsed to double.");}       
			}
			else if (peakPattern.matcher(line).matches()){
				if (peaks == null) throw new IOException("in line "+i+". Peaks given before collision energy definition.");
				String[] values = line.split("\\s+");
				if (values.length < 2) throw new IOException("in line "+i+". Line seems to contain a peak, but does not.");
				try { peaks.add(new Peak(Double.parseDouble(values[0]), Double.parseDouble(values[1]), 0)); } // energy is 0 for now. Set later!
				catch (NumberFormatException e) {
					throw new IOException("in line "+i+". String \""+values[0]+"\" or \""+values[1]+"\" could not be parsed to double");
				}
			}
			else if (line.contains("ms1")||line.contains("mslevel 1")) {
				for (line = file.readLine(); line != null && !line.isEmpty(); line = file.readLine()){
					++i;
					if (!line.contains("#")&&!line.contains("parentmass")&&!line.contains("charge")&&!line.contains("collision")&&!line.contains("tic")) {
						int offset = 0; //line.contains("#")?1:0;
						str = line.substring(offset, line.indexOf(' '));
						double mass = Double.parseDouble(str);
						str = line.substring(line.indexOf(' ') + 1, line.length());
						double intensity = Double.parseDouble(str);
						ms1peaks.add(new Peak(mass, intensity, 0));
					}					
				}
				ms1ParentMass = ms1peaks.get(0).getMass();
				System.out.println("MS1 read");
			} else if (line.equals("") || line.charAt(0)== '#') {
				// do nothing
			} else {
				System.out.println("Input file line "+i+": \""+line+"\" ignored.");
			}
			++i;
		}

		if (tic == 0.0 && param.raw) throw new IOException("in line "+i+". TIC of previous spectrum undefined, cannot calculate raw intensities.\nUse \"-raw off\" switch or correct input file.");		
		spectra.add(new Spectrum(collisionEnergy, tic, peaks));

		Collections.sort(spectra);

		// Set the energy values
		i = 0;
		for (Spectrum s : spectra){
			for (Peak p : s.getPeaks()) {
				p.setEnergy(i);
			}
			++i;
		}

		// recalibrate masses
		/*double parentMass = 0.0, minDiff = Double.POSITIVE_INFINITY;
		if (ms1ParentMass != 0.0){
			for (Spectrum s : spectra){
				for (Peak p : s.getPeaks()){
					if (Math.abs(p.getMass()-ms1ParentMass) < minDiff){
						minDiff = Math.abs(p.getMass()-ms1ParentMass);
						parentMass = p.getMass();
					}
				}
			}
			if (parentMass != 0) {
				double factor = ms1ParentMass / parentMass;
				System.out.println("Calibration factor: " + factor);
				for (Spectrum s : spectra) {
					for (Peak p : s.getPeaks()) {
						p.setMass(p.getMass() * factor);
					}
				}
			}			
		}	*/	

		double mass = (formula != null)?formula.getMass():0.0;
		if (focusedMass == 0.0) focusedMass = mass;
		Compound result = new Compound(formula, name, mode, mass, focusedMass, "", spectra, ms1peaks); //instrument is unknown.
		if (formula2 != null) result.setFormula2(formula2);
		return result;
	}

	static void scorePMD(Molecule PMD, Parameters param) throws ElementNotFoundException{
		// Score according to H/C-Ratio
		if (param.hcRatioScoring) {
			double hcRatio = ((double) PMD.getElementAbundance("H"))/PMD.getElementAbundance("C");
			PMD.logAndAddScore(MathUtils.pdf(hcRatio, param.hcAverage, param.hcDev));
		}

		if (param.heteroScoring){
			int hetero = PMD.getElementAbundance("N")+PMD.getElementAbundance("O")+PMD.getElementAbundance("P")+PMD.getElementAbundance("S");
			double heteroRatio = ((double) hetero/PMD.getElementAbundance("C"));
			PMD.logAndAddScore(MathUtils.pdf(heteroRatio, param.heteroAverage, param.heteroDev));
		}

		// Score according to DBE distribution
		if (param.DBEDependentScoring){
			PMD.logAndAddScore(MathUtils.pdf(PMD.getDBE(), param.DBEAverage, param.DBEStdDev));
		}

		/*if (PMD.calcDBE() < 0){
    	PMD.addScore(0.45-0.1*PMD.calcDBE());
    }	*/
	}

}
