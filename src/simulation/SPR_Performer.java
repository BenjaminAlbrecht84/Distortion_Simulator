package simulation;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Random;

import util.MyNode;
import util.MyTree;

public class SPR_Performer {

	public ArrayList<BitSet[]> run(MyTree modelTree, ArrayList<MyTree> trees, int numOfSPRTrees, int numOfSPRs, ArrayList<String> taxaOrdering) {

		ArrayList<BitSet[]> appliedSPRClusters = new ArrayList<BitSet[]>();
		ArrayList<BitSet> reprClusters = new ArrayList<BitSet>();
		for (int k = 0; k < numOfSPRs; k++) {

			// getting cluster information
			// int cardinality_cutOff = (int) Math.round(new Double(taxaOrdering.size()) / 10.);
			// cardinality_cutOff = cardinality_cutOff == 1 ? 2 : cardinality_cutOff;
			int cardinality_cutOff = 3;
			HashMap<BitSet, Integer> cluster2occ = new HashMap<BitSet, Integer>();
			for (MyTree t : trees) {
				for (BitSet b : t.getAllCluster(taxaOrdering)) {
					if (b.cardinality() > cardinality_cutOff) {
						cluster2occ.putIfAbsent(b, 0);
						cluster2occ.put(b, cluster2occ.get(b) + 1);
					}
				}
			}

			// choosing rSPR clusters
			ArrayList<BitSet> allClusters = new ArrayList<BitSet>(cluster2occ.keySet());
			ArrayList<Object[]> allSPRMoves = new ArrayList<Object[]>();
			for (int i = 0; i < allClusters.size(); i++) {
				BitSet b1 = allClusters.get(i);
				if (cluster2occ.get(b1) >= numOfSPRTrees) {
					for (int j = i + 1; j < allClusters.size(); j++) {
						BitSet b2 = allClusters.get(j);
						BitSet b = (BitSet) b1.clone();
						b.and(b2);
						if (cluster2occ.get(b2) >= numOfSPRTrees && b.cardinality() == 0) {
							MyNode v1 = modelTree.getNode(b1, taxaOrdering);
							MyNode v2 = modelTree.getNode(b2, taxaOrdering);
							if (v1 != null && v2 != null) {
								double d = cmpDistance(v1, v2);
								Object[] sprMove = { b1, b2, d };
								allSPRMoves.add(sprMove);
							}
						}
					}
				}
			}

			Collections.sort(allSPRMoves, new SPR_Comparator());
			Object[] spr = null;
			for (Object[] o : allSPRMoves) {
				BitSet b1 = (BitSet) o[0], b2 = (BitSet) o[1];
				if (!reprClusters.contains(b1) && !reprClusters.contains(b2)) {
					reprClusters.add(b1);
					reprClusters.add(b2);
					spr = o;
					break;
				}
			}

			if (spr == null)
				break;

			int counter = 0;
			ArrayList<MyTree> treeSet = (ArrayList<MyTree>) trees.clone();
			while (counter != numOfSPRTrees && !treeSet.isEmpty()) {
				MyTree t = pickRandomTree(treeSet);
				MyNode v1 = t.getNode((BitSet) spr[0], taxaOrdering);
				MyNode v2 = t.getNode((BitSet) spr[1], taxaOrdering);
				if (v1 != null && v2 != null) {
					System.out.println("PickedTree: " + trees.indexOf(t));
					System.out.println(t.getCluster(v1.getParent(), taxaOrdering) + "," + spr[0] + " -> " + spr[1]);
					System.out.println(t.toNewickString());
					t.pruneAndRegraft(v1, v2, taxaOrdering);
					counter++;
					System.out.println(t.toNewickString());
					BitSet[] sprClusters = { (BitSet) spr[0], (BitSet) spr[1] };
					appliedSPRClusters.add(sprClusters);
				} else
					treeSet.remove(t);
			}

		}

		return appliedSPRClusters;

	}

	private MyTree pickRandomTree(ArrayList<MyTree> trees) {
		Random rand = new Random();
		int randIndex = rand.nextInt(trees.size());
		MyTree t = trees.get(randIndex);
		return t;
	}

	private double cmpDistance(MyNode v1, MyNode v2) {

		ArrayList<MyNode> vec1 = new ArrayList<MyNode>();
		cmpAncestorsRec(v1, vec1);
		ArrayList<MyNode> vec2 = new ArrayList<MyNode>();
		cmpAncestorsRec(v2, vec2);
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

	private void cmpAncestorsRec(MyNode v, ArrayList<MyNode> vec) {
		MyNode p = v.getParent();
		if (p != null) {
			vec.add(p);
			cmpAncestorsRec(p, vec);
		}
	}

	public class SPR_Comparator implements Comparator<Object[]> {

		@Override
		public int compare(Object[] o1, Object[] o2) {
			double d1 = (double) o1[2];
			double d2 = (double) o2[2];
			if (d1 > d2)
				return -1;
			else if (d2 > d1)
				return 1;
			return 0;
		}

	}

}
