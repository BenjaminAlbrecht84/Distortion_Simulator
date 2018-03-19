package method;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Random;

import algorithm.linear.tree.LinearMajorityTree;
import util.MyNewickParser;
import util.MyNode;
import util.MyTree;
import util.graph.MyPhyloTree;
import util.newick.NewickParser;

public class Distortion_Manager {

	public void run(MyTree speciesTree, ArrayList<MyTree> geneTrees, ArrayList<String> taxaOrdering, boolean addUniqueTaxa, int startIndex) {

		// calculating consensus tree
		MyTree consensusTree = speciesTree;
		if (consensusTree == null) {
			ArrayList<MyPhyloTree> phyloTrees = new ArrayList<MyPhyloTree>();
			for (MyTree t : geneTrees)
				phyloTrees.add(new MyPhyloTree(new NewickParser().run(t.toNewickString())));
			MyPhyloTree consensusPyloTree = new LinearMajorityTree(taxaOrdering).run(phyloTrees, 0.51);
			consensusTree = new MyNewickParser().run(consensusPyloTree.toBracketString());
		}

		// adding unique taxa
		if (addUniqueTaxa) {
			for (MyTree t : geneTrees) {
				for (int i = 0; i < 2; i++)
					addUniqueTaxa(t, taxaOrdering);
			}
		}

		// reporting output
		System.out.println(">GeneTrees: ");
		System.out.println(consensusTree.toNewickString());
		for (MyTree t : geneTrees)
			System.out.println(t.toNewickString());

		// analyzing trees
		// ArrayList<Object[]> conflictingTrees = new Tree_Analyzer().run(consensusTree, geneTrees, taxaOrdering);
		ArrayList<ArrayList<Object[]>> conflictingTrees = new Tree_Analyzer().run(consensusTree, geneTrees, taxaOrdering, startIndex);

		System.out.println("ConflictingTrees: ");
		int i1 = 0, i2 = 0;
		double max = 0, sum = 0;
		int counter = 0;
		for (ArrayList<Object[]> list : conflictingTrees) {

			for (Object[] o : list) {
				if ((double) o[2] > max) {
					counter++;
					sum += (double) o[2];
					max = (double) o[2];
					i1 = conflictingTrees.indexOf(list);
					i2 = list.indexOf(o);
				}

				System.out.println("> " + o[2] + " --- rSPR-Moves: " + o[1]);
				System.out.print("HGT-Sets: ");
				for (BitSet s : (ArrayList<BitSet>) o[3])
					System.out.print(s + " ");
				System.out.println("\n" + consensusTree.toNewickString());
				System.out.println(((MyTree) o[0]).toNewickString());
				
			}

		}

		System.out.println("-------- ");
		Object[] o = conflictingTrees.get(i1).get(i2);
		System.out.println("> Dist: " + o[2] + " (" + Math.round((sum-max) / (counter-1)) + "), rSPR-Moves: " + o[1]);
		System.out.print("HGT-Sets: ");
		for (BitSet s : (ArrayList<BitSet>) o[3])
			System.out.print(s + " ");
		System.out.println("\n" + consensusTree.toNewickString());
		System.out.println(((MyTree) o[0]).toNewickString());

	}

	private void addUniqueTaxa(MyTree t, ArrayList<String> taxaOrdering) {
		int id = 1;
		while (taxaOrdering.contains(String.valueOf(id)))
			id++;
		taxaOrdering.add(String.valueOf(id));
		MyNode l = new MyNode(String.valueOf(id));
		ArrayList<MyNode> nodes = t.getNodes();
		Random rand = new Random();
		MyNode r = null;
		while (r == null) {
			MyNode v = nodes.get(rand.nextInt(nodes.size()));
			if (!v.isLeaf()) {
				v.addChildren(l);
				break;
			}
		}
	}

}
