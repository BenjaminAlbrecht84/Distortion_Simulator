package method;

import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;

import util.MyNewickParser;
import util.MyNode;
import util.MyTree;
import util.Splits_Calculator;
import util.TreeRooter;

public class Tree_Analyzer {

	private MyDistortion_Scorer distScorer;

	public ArrayList<ArrayList<Object[]>> run(MyTree speciesTree, ArrayList<MyTree> geneTrees, ArrayList<String> taxaOrdering, int startIndex) {

		ArrayList<ArrayList<Object[]>> conflictingTreeSet = new ArrayList<ArrayList<Object[]>>();
		for (int i = startIndex; i < geneTrees.size(); i++) {

			MyTree geneTree = geneTrees.get(i);

			System.out.println(geneTrees.indexOf(geneTree) + " / " + geneTrees.size());
			System.out.println(">\n" + speciesTree.toNewickString());
			System.out.println(geneTree.toNewickString());

			MyTree speciesTreeRestricted = speciesTree.copy();
			MyTree geneTreeRestricted = geneTree.copy();
			ArrayList<Object[]> replacementInfo = adaptTaxa(speciesTreeRestricted, geneTreeRestricted, taxaOrdering);
			ArrayList<String> taxaOrderingRestricted = speciesTreeRestricted.getTaxa();
			Collections.sort(taxaOrderingRestricted, new TaxaComparator(taxaOrdering));

			distScorer = new MyDistortion_Scorer(speciesTreeRestricted, taxaOrderingRestricted);

			MyTree resolvedGeneTree = new MyNewickParser().run(geneTreeRestricted.toNewickString());

			// resolving distortion in gene tree

			HashMap<MyNode, Integer> node2distortion = initDistScore(speciesTreeRestricted, resolvedGeneTree, taxaOrderingRestricted);
			ArrayList<Object[]> distortionHistory = new ArrayList<Object[]>();
			ArrayList<MyNode> endPoints;
			while (!(endPoints = cmpEndPoints(node2distortion, resolvedGeneTree, taxaOrderingRestricted)).isEmpty()) {
				resolveDistortion(speciesTreeRestricted, resolvedGeneTree, endPoints, node2distortion, taxaOrderingRestricted, distortionHistory);
				node2distortion = initDistScore(speciesTreeRestricted, resolvedGeneTree, taxaOrderingRestricted);
			}

			ArrayList<Object[]> conflictingTrees = Distortion_Adder.run(speciesTreeRestricted, geneTreeRestricted, resolvedGeneTree,
					taxaOrderingRestricted, taxaOrdering, replacementInfo, distortionHistory);
			conflictingTreeSet.add(conflictingTrees);

		}

		return conflictingTreeSet;

	}

	private ArrayList<Object[]> adaptTaxa(MyTree speciesTree, MyTree geneTree, ArrayList<String> taxaOrdering) {

		ArrayList<String> speciesTreeTaxa = speciesTree.getTaxa();
		ArrayList<String> geneTreeTaxa = geneTree.getTaxa();

		ArrayList<MyNode> uniqueSpeciesLeaves = getUniqueLeaves(speciesTree, geneTreeTaxa);
		ArrayList<MyNode> uniqueGeneLeaves = getUniqueLeaves(geneTree, speciesTreeTaxa);

		removeUniqueLeaves(uniqueSpeciesLeaves, speciesTree, taxaOrdering);
		ArrayList<Object[]> replacementInfo = removeUniqueLeaves(uniqueGeneLeaves, geneTree, taxaOrdering);

		return replacementInfo;
	}

	private void reRootSpeciesTree(MyTree speciesTree, MyTree geneTree, ArrayList<String> taxaOrdering) {
		TreeRooter.run(speciesTree, geneTree, taxaOrdering);
		distScorer = new MyDistortion_Scorer(speciesTree, taxaOrdering);
	}

	private ArrayList<Object[]> removeUniqueLeaves(ArrayList<MyNode> uniqueLeaves, MyTree t, ArrayList<String> taxaOrdering) {
		ArrayList<Object[]> replacementInfo = new ArrayList<Object[]>();
		for (MyNode l : uniqueLeaves) {
			String taxon = l.getId();
			BitSet cluster = t.getCluster(l.getParent(), taxaOrdering);
			cluster.set(taxaOrdering.indexOf(taxon), false);
			boolean addAbove = l.getParent().getChildren().size() == 2;
			MyNode p = l.getParent();
			t.removeNode(l, false, taxaOrdering);
			double[] branchLengths = { p.getBranchLength(), p.getChildren().get(0).getBranchLength() };
			t.supressNode(p, taxaOrdering);
			Object[] o = { taxon, cluster, addAbove, branchLengths };
			replacementInfo.add(o);
		}
		return replacementInfo;
	}

	private ArrayList<MyNode> getUniqueLeaves(MyTree t, ArrayList<String> taxaSet) {
		ArrayList<MyNode> uniqueLeaves = new ArrayList<MyNode>();
		for (MyNode l : t.getLeafNodes()) {
			if (!taxaSet.contains(l.getId()))
				uniqueLeaves.add(l);
		}
		return uniqueLeaves;
	}

	private HashMap<MyNode, Integer> initDistScore(MyTree t, MyTree geneTree, ArrayList<String> taxaOrdering) {
		ArrayList<BitSet[]> splits = Splits_Calculator.run(geneTree, taxaOrdering);
		HashMap<MyNode, Integer> node2distortion = new HashMap<MyNode, Integer>();
		for (BitSet[] split : splits) {

			int distortion = 0;
			try {
				distortion = distScorer.run(split[0], split[1]);
			} catch (IOException e) {
				distortion = -1;
				e.printStackTrace();
			}

			if (distortion != 0 && distortion != -1) {
				MyNode v = geneTree.getNode(split[0], taxaOrdering);
				node2distortion.put(v, distortion);
			}

		}
		return node2distortion;
	}

	private void resolveDistortion(MyTree speciesTree, MyTree geneTree, ArrayList<MyNode> endPoints, HashMap<MyNode, Integer> node2distortion,
			ArrayList<String> taxaOrdering, ArrayList<Object[]> distortionHistory) {

		boolean distortionDecreased = false;
		Collections.sort(endPoints, new EndpointComaprator(geneTree, taxaOrdering));
		for (MyNode v : endPoints) {

			// System.out.println(">" + geneTree.getCluster(v, taxaOrdering));
			// System.out.println(geneTree.toNewickString());
			resolveTreeRec(speciesTree, geneTree, taxaOrdering, node2distortion);
			// System.out.println(geneTree.toNewickString());

			Object[] rsprSet = selectBestRSPRMove(v, speciesTree, geneTree, taxaOrdering, node2distortion);
			// System.out.println(geneTree.toNewickString());

			if (rsprSet != null) {

				// System.out.println(rsprSet[2]);

				// removing distortion by applying selected rSPR moves
				MyNode v1 = geneTree.getNode((BitSet) rsprSet[0], taxaOrdering);
				MyNode v2 = geneTree.getNode((BitSet) rsprSet[1], taxaOrdering);

				// System.out.println(rsprSet[0] + " -> " + rsprSet[1] + " : " + rsprSet[2]);
				// System.out.println(geneTree.toNewickString());

				BitSet v1Set = (BitSet) rsprSet[0];
				MyNode p1 = v1.getParent();
				geneTree.removeNode(v1, false, taxaOrdering);
				v2.addChildren(v1);
				BitSet p1Set = geneTree.getCluster(p1, taxaOrdering);

				double[] branchLengths = new double[2];
				if (p1.getChildren().size() == 1) {
					branchLengths[0] = p1.getBranchLength();
					branchLengths[1] = p1.getChildren().get(0).getBranchLength();
				}

				boolean addAbove = geneTree.supressNode(p1, taxaOrdering);

				Object[] distHistory = { v1Set, p1Set, addAbove, branchLengths };
				distortionHistory.add(distHistory);

				distortionDecreased = true;
				break;

			}

		}

		if (!distortionDecreased)
			reRootSpeciesTree(speciesTree, geneTree, taxaOrdering);

		// System.out.println(geneTree.toNewickString());

	}

	private void resolveTreeRec(MyTree t, MyTree geneTree, ArrayList<String> taxaOrdering, HashMap<MyNode, Integer> node2distortion) {
		ArrayList<MyNode> nodes = geneTree.getNodes();
		boolean isResolved = false;
		for (MyNode v : nodes) {
			if (resolveNode(t, geneTree, v, taxaOrdering, node2distortion))
				isResolved = true;
		}
		if (isResolved)
			resolveTreeRec(t, geneTree, taxaOrdering, node2distortion);
	}

	private boolean resolveNode(MyTree t, MyTree geneTree, MyNode v, ArrayList<String> taxaOrdering, HashMap<MyNode, Integer> node2distortion) {

		if (v.getChildren().size() < 3)
			return false;

		ArrayList<BitSet> childSets = new ArrayList<BitSet>();
		for (MyNode w : v.getChildren())
			childSets.add(geneTree.getCluster(w, taxaOrdering));

		ArrayList<BitSet> toResolve = new ArrayList<BitSet>();
		BitSet consideredCluster = new BitSet(taxaOrdering.size());
		for (MyNode c : v.getChildren()) {
			if (node2distortion.containsKey(c))
				continue;
			BitSet cCluster = geneTree.getCluster(c, taxaOrdering);
			if (containmentCheck(consideredCluster, cCluster))
				continue;
			MyNode w = t.getNode(cCluster, taxaOrdering);
			if (w != null) {
				BitSet wCluster = t.getCluster(w, taxaOrdering);
				BitSet cluster = wCluster;
				while (w.getParent() != null) {
					w = w.getParent();
					BitSet b = (BitSet) t.getCluster(w, taxaOrdering).clone();
					b.xor(cluster);
					if (!childSets.contains(b))
						break;
					cluster = (BitSet) t.getCluster(w, taxaOrdering).clone();
				}
				consideredCluster.or(cluster);
				toResolve.add(cluster);
			}
		}

		boolean isResolved = false;
		for (BitSet b : toResolve) {
			if (establishChildCluster(geneTree, v, b, taxaOrdering))
				isResolved = true;
		}

		return isResolved;

	}

	private static boolean establishChildCluster(MyTree geneTree, MyNode v, BitSet cluster, ArrayList<String> taxaOrdering) {

		ArrayList<MyNode> toRefine = new ArrayList<MyNode>();
		for (MyNode w : v.getChildren()) {
			BitSet c = geneTree.getCluster(w, taxaOrdering);
			BitSet b = (BitSet) cluster.clone();
			b.or(c);
			if (b.equals(cluster))
				toRefine.add(w);
		}
		boolean isRefined = false;
		if (toRefine.size() < v.getChildren().size() && toRefine.size() > 1) {

			// System.out.println(cluster);
			// System.out.println(geneTree.toNewickString());

			MyNode x = new MyNode();
			v.addChildren(x);
			for (MyNode y : toRefine) {
				v.removeChild(y);
				x.addChildren(y);
			}
			isRefined = true;

			// System.out.println(geneTree.toNewickString());
		}

		return isRefined;

	}

	private Object[] selectBestRSPRMove(MyNode endPoint, MyTree t, MyTree geneTree, ArrayList<String> taxaOrdering,
			HashMap<MyNode, Integer> node2distortion) {

		// System.out.println();

		try {

			int dist1 = 0;
			for (MyNode v : geneTree.getNodes())
				dist1 += distScorer.run(geneTree.getCluster(v, taxaOrdering));

			int min = dist1;
			Object[] result = null;
			for (int i = 0; i < geneTree.getNodes().size(); i++) {

				MyNode w = geneTree.getNodes().get(i);
				MyNode[] nodes = { endPoint, w };

				for (MyNode x : nodes) {

					MyNode y = nodes[0].equals(x) ? nodes[1] : nodes[0];

					ArrayList<MyNode> children = new ArrayList<MyNode>();
					for (MyNode c : x.getChildren()) {
						if (!node2distortion.containsKey(c))
							children.add(c);
					}
					// if (x.getParent() != null && x.getParent().equals(geneTree.getRoot())) {
					// for (MyNode c : geneTree.getRoot().getChildren()) {
					// if (!w.equals(x) && !node2distortion.containsKey(c))
					// children.add(c);
					// }
					// }

					for (MyNode c : children) {

						BitSet cCluster = geneTree.getCluster(c, taxaOrdering);
						BitSet yCluster = geneTree.getCluster(y, taxaOrdering);
						if (containmentCheck(cCluster, yCluster) || c.getParent().equals(y))
							continue;

						MyTree geneTreeCopy = geneTree.copy();
						MyNode cCopy = geneTreeCopy.getNodes().get(geneTree.getNodes().indexOf(c));
						boolean isRootChild = cCopy.getParent().equals(geneTreeCopy.getRoot());

						MyNode yCopy = geneTreeCopy.getNodes().get(geneTree.getNodes().indexOf(y));
						geneTreeCopy.removeNode(cCopy, true, taxaOrdering);
						if (yCopy.isLeaf())
							geneTreeCopy.attachNodeAbove(cCopy, yCopy, null);
						else
							yCopy.addChildren(cCopy);

						int dist2 = 0;
						for (MyNode v : geneTreeCopy.getNodes())
							dist2 += distScorer.run(geneTreeCopy.getCluster(v, taxaOrdering));

						// System.out.println(cCluster + " " + yCluster + " " + dist2 + " > " + dist1);
						if (dist2 < min) { // || (dist2 == min && isRootChild)) {
							min = dist2;
							result = new Object[3];
							result[0] = cCluster;
							result[1] = yCluster;
							result[2] = dist2;
						}

					}

				}

			}

			// if (result == null) {
			// System.out.println("ERROR - No decrease in distortion: " + result[2] + " >= " + dist1);
			// System.out.println(t.toNewickString());
			// System.out.println(geneTree.toNewickString());
			// System.out.println(geneTree.getCluster(endPoint, taxaOrdering) + ": " + result[0] + " -> " + result[1]);
			// // System.exit(0);
			// }

			return result;

		} catch (Exception e) {
			e.printStackTrace();
		}

		return null;

	}

	private boolean containmentCheck(BitSet b1, BitSet b2) {
		BitSet b = (BitSet) b1.clone();
		b.or(b2);
		return b.equals(b1);
	}

	private ArrayList<MyNode> cmpEndPoints(HashMap<MyNode, Integer> node2distortion, MyTree geneTree, ArrayList<String> taxaOrdering) {
		ArrayList<MyNode> endPoints = new ArrayList<MyNode>();
		for (MyNode v : node2distortion.keySet()) {
			int counter = 0;
			for (MyNode w : v.getChildren()) {
				if (!node2distortion.containsKey(w))
					counter++;
			}
			if (counter != 0)
				endPoints.add(v);
		}
		return endPoints;
	}

	public class TaxaComparator implements Comparator<String> {

		private ArrayList<String> taxaOrdering;

		public TaxaComparator(ArrayList<String> taxaOrdering) {
			this.taxaOrdering = taxaOrdering;
		}

		@Override
		public int compare(String s1, String s2) {
			int p1 = taxaOrdering.indexOf(s1);
			int p2 = taxaOrdering.indexOf(s2);
			if (p1 < p2)
				return -1;
			else if (p1 > p2)
				return 1;
			return 0;
		}

	}

	public class EndpointComaprator implements Comparator<MyNode> {

		private MyTree t;
		private ArrayList<String> taxaOrdering;

		public EndpointComaprator(MyTree t, ArrayList<String> taxaOrdering) {
			this.t = t;
			this.taxaOrdering = taxaOrdering;
		}

		@Override
		public int compare(MyNode v1, MyNode v2) {
			BitSet b1 = t.getCluster(v1, taxaOrdering);
			BitSet b2 = t.getCluster(v2, taxaOrdering);
			int i1 = b1.nextSetBit(0);
			int i2 = b2.nextSetBit(0);
			while (i1 == i2) {
				i1 = b1.nextSetBit(i1 + 1);
				i2 = b2.nextSetBit(i2 + 1);
			}
			return Integer.compare(i1, i2);
		}

	}

}
