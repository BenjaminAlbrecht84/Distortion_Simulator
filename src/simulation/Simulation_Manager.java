package simulation;

import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Random;

import util.MyNode;
import util.MyTree;
import util.RandomBinaryTree_Generator;

public class Simulation_Manager {

	private final static char[] SIGMA = { 'A', 'G', 'T', 'C' };
	private double maxPE = 0., minPE = Double.POSITIVE_INFINITY, avgPE = 0, peCounter = 0;

	public Object[] run(int numberOfTaxa, int maxBranchLength, int seqLength, int numOfTrees, int numOfSPRTrees, int numOfSPRs, double u) throws IOException {

		// creating random tree
		System.out.println(">Step 1 - Initializing Jukes Cantor Model");
		Object[] res = RandomBinaryTree_Generator.run(numberOfTaxa);
		MyTree speciesTree = (MyTree) res[0];
		ArrayList<String> taxaOrdering = (ArrayList<String>) res[1];
		speciesTree.assingBranchLength(maxBranchLength, numberOfTaxa);
		String startSequence = createRandomSequence(seqLength);
		System.out.println(speciesTree.toNewickString());
		System.out.println("Avg-EdgeLength: " + speciesTree.getAvgEdgeLength());

		// computing sets of evolved sequences
		System.out.println(">Step 2 - Evolving sequences along the model tree");
		ArrayList<HashMap<MyNode, String>> setOfEvolvedSeqs = evolveSequence(speciesTree, startSequence, numOfTrees, u);
		System.out.println("MaxPE: " + maxPE + ", minPE: " + minPE + " AvgPE: " + avgPE);
		System.out.println(setOfEvolvedSeqs.size() + "x" + setOfEvolvedSeqs.get(0).size() + " sequences generated.");

		// inferring trees based on evolved sequences
		System.out.println(">Step 3 - Inferring trees from evolved sequences");
		ArrayList<MyTree> geneTrees = new RAxML_TreeGenerator().run(setOfEvolvedSeqs);
		System.out.println(geneTrees.size() + " trees inferred.");
		for (MyTree t : geneTrees)
			System.out.println(t.toNewickString());

		// applying rSPR-moves
		System.out.println(">Step 4 - Applying rSPR moves");
		ArrayList<BitSet[]> sprMoves = new SPR_Performer().run(speciesTree, geneTrees, numOfSPRTrees, numOfSPRs, taxaOrdering);
		for (Object[] spr : sprMoves)
			System.out.println("SPR-Clusters: " + spr[0] + " " + spr[1]);
		System.out.println(numOfSPRs + "x" + numOfSPRTrees + " tree(s) modified.");

		// reporting output
		System.out.println(">Output ");
		for (MyTree t : geneTrees)
			System.out.println(t.toNewickString());

		Object[] result = { speciesTree, geneTrees, taxaOrdering };
		return result;

	}



	private ArrayList<HashMap<MyNode, String>> evolveSequence(MyTree modelTree, String startSequence, int numOfTrees, double u) {
		ArrayList<HashMap<MyNode, String>> setOfevolvedSeqs = new ArrayList<HashMap<MyNode, String>>();
		for (int i = 0; i < numOfTrees; i++) {
			HashMap<MyNode, String> evolvedSeqs = new HashMap<MyNode, String>();
			evolveSeqRec(modelTree.getRoot(), startSequence, evolvedSeqs, u);
			setOfevolvedSeqs.add(evolvedSeqs);
		}
		avgPE /= peCounter;

		return setOfevolvedSeqs;
	}

	private void evolveSeqRec(MyNode v, String sequence, HashMap<MyNode, String> evolvedSeqs, double u) {
		if (v.isLeaf())
			evolvedSeqs.put(v, sequence);
		else {
			for (MyNode c : v.getChildren()) {
				StringBuilder evolvedSeq = new StringBuilder();
				for (int i = 0; i < sequence.length(); i++) {
					double t = c.getBranchLength();
					double p = 1. - Math.pow(Math.E, (-4. / 3.) * u * t);
					if (Math.random() < p)
						evolvedSeq.append(randomChar());
					else
						evolvedSeq.append(sequence.charAt(i));
					maxPE = p > maxPE ? p : maxPE;
					minPE = p < minPE ? p : minPE;
					avgPE += p;
					peCounter++;
				}
				evolveSeqRec(c, evolvedSeq.toString(), evolvedSeqs, u);
			}
		}
	}

	private String createRandomSequence(int seqLength) {
		StringBuilder s = new StringBuilder();
		for (int i = 0; i < seqLength; i++)
			s.append(randomChar());
		return s.toString();
	}

	private char randomChar() {
		double d = Math.random() * (SIGMA.length - 1);
		int index = (int) Math.round(d);
		return SIGMA[index];
	}

}
