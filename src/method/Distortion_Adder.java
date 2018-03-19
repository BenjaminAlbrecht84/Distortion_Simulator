package method;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;

import util.MyNode;
import util.MyTree;
import util.TreeRooter;

public class Distortion_Adder {

	public static ArrayList<Object[]> run(MyTree speciesTree, MyTree geneTree, MyTree resolvedGeneTree, ArrayList<String> taxaOrderingRestricted,
			ArrayList<String> taxaOrdering, ArrayList<Object[]> replacementInfo, ArrayList<Object[]> distortionHistory) {

//		System.out.println("---------");
//		System.out.println(speciesTree.toNewickString());
//		System.out.println(geneTree.toNewickString());
//		System.out.println(resolvedGeneTree.toNewickString());

		MyTree resolvedGeneTreeCopy = resolvedGeneTree.copy();

		// creating distortion tree with distortion links
		MyTree distortionTree = speciesTree.copy();
		TreeRooter.run(distortionTree, resolvedGeneTreeCopy, taxaOrderingRestricted);

		// System.out.println(distortionTree.toNewickString());

		addDistortionLinks(distortionTree, resolvedGeneTreeCopy, taxaOrderingRestricted, distortionHistory);

		// System.out.println(distortionTree.toNewickString());

		addUniqueTaxa(distortionTree, replacementInfo, taxaOrdering);

		// System.out.println(distortionTree.toNewickString());

		// compute distances
		cmpTreeDistances(distortionTree, taxaOrdering);
		ArrayList<Double> allDistances = new ArrayList<Double>();
		for (MyNode v : distortionTree.getNodes()) {
			if (v.getDistortionOperation() != null) {
				allDistances.add((Double) v.getDistortionOperation()[2]);
			}
		}
		allDistances.add(0.);
		Collections.sort(allDistances);

		// generating filtered distortion trees
		ArrayList<Object[]> conflictingTrees = new ArrayList<Object[]>();
		// for (double distTreshold : allDistances) {
		for (int i = 0; i < allDistances.size(); i++) {
			double distTreshold = allDistances.get(i);
			MyTree distortionTreeCopy = distortionTree.copy();
			Object[] filteringResult = filterDistortion(distortionTreeCopy, distTreshold, taxaOrdering);
			if ((int) filteringResult[1] > 0) {
				MyTree conflictingTree = (MyTree) filteringResult[0];
				if (conflictingTree.getLastRoot() != null)
					TreeRooter.run(conflictingTree, conflictingTree.getLastRoot(), taxaOrdering);
				Object[] distResult = { conflictingTree, filteringResult[1], allDistances.get(i + 1), filteringResult[2] };
				conflictingTrees.add(distResult);
			} else
				break;
		}

		return conflictingTrees;

	}

	private static void addUniqueTaxa(MyTree t, ArrayList<Object[]> replacementInfo, ArrayList<String> taxaOrdering) {
		for (int i = replacementInfo.size() - 1; i >= 0; i--) {
			Object[] o = replacementInfo.get(i);
			MyNode l = new MyNode((String) o[0]);
			MyNode v = t.getLCA((BitSet) o[1], taxaOrdering);
			double[] brachLengths = cmpBranchLengths(v, (double[]) o[3]);
			if (!(boolean) o[2])
				v.addChildren(l);
			else
				t.attachNodeAbove(l, v, brachLengths);
		}
	}

	private static Object[] filterDistortion(MyTree distortionTree, double distTreshold, ArrayList<String> taxaOrdering) {

		MyTree conflictingTree = distortionTree.copy();
		cmpTreeDistances(conflictingTree, taxaOrdering);
		ArrayList<MyNode> distOrdering = cmpDistOrdering(conflictingTree);
		int counter = distOrdering.size();

		ArrayList<MyNode> abandonedNodes = new ArrayList<MyNode>();
		ArrayList<MyNode> filteredNodes = new ArrayList<MyNode>();
		for (MyNode v : distOrdering) {
			if ((Double) v.getDistortionOperation()[2] <= distTreshold && distortionCheck(v)) {
				abandonedNodes.add(v.getParent());
				conflictingTree.removeNode(v, false, taxaOrdering);
				((MyNode) v.getDistortionOperation()[0]).addChildren(v);
				counter--;
			} else
				filteredNodes.add(v);
		}

		ArrayList<BitSet> filteredSets = new ArrayList<BitSet>();
		for (MyNode v : filteredNodes)
			filteredSets.add(conflictingTree.getCluster(v, taxaOrdering));

		cleanUpRec(conflictingTree, taxaOrdering);

		Object[] result = { conflictingTree, counter, filteredSets };
		return result;
	}

	private static boolean distortionCheck(MyNode v) {
		MyNode w = (MyNode) v.getDistortionOperation()[0];
		MyNode p = w.getParent();
		while (p != null && !p.equals(v))
			p = p.getParent();
		if (p != null && p.equals(v))
			return false;
		return true;
	}

	private static ArrayList<MyNode> cmpDistOrdering(MyTree t) {
		ArrayList<MyNode> distortionOperations = new ArrayList<MyNode>();
		for (MyNode v : t.getNodes()) {
			if (v.getDistortionOperation() != null)
				distortionOperations.add(v);
		}
		Collections.sort(distortionOperations, new Distortion_Adder.DistortionComparator());
		return distortionOperations;
	}

	static class DistortionComparator implements Comparator<MyNode> {

		@Override
		public int compare(MyNode v1, MyNode v2) {
			int i1 = (int) v1.getDistortionOperation()[1];
			int i2 = (int) v2.getDistortionOperation()[1];
			return Integer.compare(i1, i2);
		}

	}

	private static void cleanUpRec(MyTree t, ArrayList<String> taxaOrdering) {
		ArrayList<MyNode> isolatedNodes = new ArrayList<MyNode>();
		ArrayList<MyNode> nonLeafNodes = new ArrayList<MyNode>();
		for (MyNode v : t.getNodes()) {
			if (v.getChildren().size() == 1)
				isolatedNodes.add(v);
			if (v.isLeaf() && v.getId().isEmpty())
				nonLeafNodes.add(v);
		}
		for (MyNode v : isolatedNodes) {
			if (v.equals(t.getLastRoot()) && v.getChildren().size() == 1) {
				t.setLastRoot(v.getChildren().get(0));
			}
			t.supressNode(v, taxaOrdering);
		}
		for (MyNode v : nonLeafNodes) {
			if (v.equals(t.getLastRoot())) {
				t.setLastRoot(v.getParent());
			}
			t.removeNode(v, false, taxaOrdering);
		}
		if (!isolatedNodes.isEmpty() || !nonLeafNodes.isEmpty())
			cleanUpRec(t, taxaOrdering);
	}

	private static void cmpTreeDistances(MyTree t, ArrayList<String> taxaOrdering) {
		for (MyNode v : t.getNodes()) {
			if (v.getDistortionOperation() != null) {
				MyNode v1 = v.getParent();
				MyNode v2 = (MyNode) v.getDistortionOperation()[0];
				double dist = cmpDistance(v1, v2, t, taxaOrdering) + 1;
				v.getDistortionOperation()[2] = dist;
			}
		}
	}

	private static double cmpDistance(MyNode v1, MyNode v2, MyTree conflictingTree, ArrayList<String> taxaOrdering) {

		ArrayList<MyNode> vec1 = new ArrayList<MyNode>();
		cmpAncestorsRec(v1, vec1, conflictingTree, taxaOrdering);
		ArrayList<MyNode> vec2 = new ArrayList<MyNode>();
		cmpAncestorsRec(v2, vec2, conflictingTree, taxaOrdering);
		MyNode lca = null;
		for (MyNode v : vec1) {
			if (vec2.contains(v)) {
				lca = v;
				break;
			}
		}

		double dist = 0;
		for (MyNode v : vec1) {
			if (v.equals(lca))
				break;
			dist += v.getBranchLength();
		}
		for (MyNode v : vec2) {
			if (v.equals(lca))
				break;
			dist += v.getBranchLength();
		}

		return dist;
	}

	private static void cmpAncestorsRec(MyNode v, ArrayList<MyNode> vec, MyTree conflictingTree, ArrayList<String> taxaOrdering) {
		if (v != null) {
			vec.add(v);
			MyNode p = v.getDistortionOperation() != null ? (MyNode) v.getDistortionOperation()[0] : v.getParent();
			cmpAncestorsRec(p, vec, conflictingTree, taxaOrdering);
		}
	}

	private static void addDistortionLinks(MyTree conflictingTree, MyTree resolvedGeneTree, ArrayList<String> taxaOrdering,
			ArrayList<Object[]> distortionHistory) {

		// adding distortion simultaneously to both conflicting tree and resolved gene tree
		// undoing all applied rSPR-moves by creating distortion links
		for (int i = distortionHistory.size() - 1; i >= 0; i--) {

			// adding distortion to resolved gene tree
			// -------------------------------------------

			Object[] o = distortionHistory.get(i);
			BitSet[] rSPR = { (BitSet) o[0], (BitSet) o[1] };
			boolean addAbove = (boolean) o[2];

			// adding distortion to resolved gene tree
			MyNode v1 = resolvedGeneTree.getNode(rSPR[0], taxaOrdering);
			MyNode v2 = resolvedGeneTree.getNode(rSPR[1], taxaOrdering);
			BitSet v1Set = resolvedGeneTree.getCluster(v1, taxaOrdering);
			BitSet v2Set = resolvedGeneTree.getCluster(v2, taxaOrdering);
			if (!v1.getParent().equals(v2) || v2.getChildren().size() > 2) {
				resolvedGeneTree.removeNode(v1, true, taxaOrdering);
				if (!addAbove)
					v2.addChildren(v1);
				else
					resolvedGeneTree.attachNodeAbove(v1, v2, null);
			}

			// -------------------------------------------

			// adding distortion to conflicting gene tree
			// -------------------------------------------

			resolveTree(conflictingTree, rSPR, taxaOrdering);

			v1 = conflictingTree.getNode(v1Set, taxaOrdering);
			v2 = conflictingTree.getNode(v2Set, taxaOrdering);

			if (!v1.getParent().equals(v2) || v2.getChildren().size() > 2) {

				MyNode p1 = v1.getParent();
				conflictingTree.removeNode(v1, false, taxaOrdering);
				if (!addAbove && !v2.isLeaf()) {
					v2.addChildren(v1);
				} else {
					double[] branchLengths = cmpBranchLengths(v2, (double[]) o[3]);
					conflictingTree.attachNodeAbove(v1, v2, branchLengths);
				}
				if (v1.getDistortionOperation() == null) {
					Object[] distortionOperation = { p1, i, null };
					v1.setDistortionOperation(distortionOperation);
				}
				if (p1.isLeaf() && p1.equals(conflictingTree.getLastRoot()))
					conflictingTree.setLastRoot(v1.getParent());

			}

			// -------------------------------------------

		}

	}

	private static double[] cmpBranchLengths(MyNode v, double[] branchLengths) {
		double sum = branchLengths[0] + branchLengths[1];
		double l1, l2;
		if (sum == 0) {
			l1 = 0;
			l2 = 0;
		} else {
			double p1 = branchLengths[0] / sum;
			double p2 = branchLengths[1] / sum;
			l1 = p1 * v.getBranchLength();
			l2 = p2 * v.getBranchLength();
		}
		double[] result = { l1, l2 };
		return result;
	}

	private static void resolveTree(MyTree t, BitSet[] clusterSet, ArrayList<String> taxaOrdering) {
		for (BitSet c : clusterSet)
			resolveTree(t, c, taxaOrdering);
	}

	private static void resolveTree(MyTree t, BitSet cluster, ArrayList<String> taxaOrdering) {

		// resolving non-binary node
		MyNode lca = t.getLCA(cluster, taxaOrdering);

		ArrayList<MyNode> toRefine = new ArrayList<MyNode>();
		for (MyNode w : lca.getChildren()) {
			BitSet c = t.getCluster(w, taxaOrdering);
			BitSet b = (BitSet) cluster.clone();
			b.or(c);
			if (b.equals(cluster))
				toRefine.add(w);
		}
		if (toRefine.size() < lca.getChildren().size() && toRefine.size() > 1) {
			MyNode x = new MyNode();
			lca.addChildren(x);
			for (MyNode y : toRefine) {
				lca.removeChild(y);
				x.addChildren(y);
			}
		}

		// re-rooting tree
		if (lca.equals(t.getRoot()) && t.getNode(cluster, taxaOrdering) == null)
			TreeRooter.run(t, cluster, taxaOrdering);

	}

}
