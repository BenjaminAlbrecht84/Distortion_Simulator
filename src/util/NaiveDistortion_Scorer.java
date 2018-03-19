package util;

import java.util.ArrayList;
import java.util.BitSet;

public class NaiveDistortion_Scorer {

	private int distortion = -1;

	public int run(MyTree t, ArrayList<String> taxaOrdering, BitSet A, BitSet B) {
		for (int k = 0; k < 5; k++) {
			cmpDistortionScoreRec(t, taxaOrdering, A, B, 0, k);
			if (distortion != -1)
				break;
		}
		return distortion;
	}

	private void cmpDistortionScoreRec(MyTree t, ArrayList<String> taxaOrdering, BitSet A, BitSet B, int score, int k) {

		if (checkTree(t, taxaOrdering, A, B))
			distortion = score;
		else if (score < k && distortion == -1) {
			ArrayList<BitSet> clusterSet = t.getAllCluster(taxaOrdering);
			for (int i = 0; i < clusterSet.size(); i++) {
				BitSet b1 = clusterSet.get(i);
				if (b1.cardinality() < taxaOrdering.size()) {
					for (int j = 0; j < clusterSet.size(); j++) {
						BitSet b2 = clusterSet.get(j);
						if (!containmentCheck(b1, b2)) {
							MyTree tCopy = new MyNewickParser().run(t.toNewickString());
							MyNode v1 = tCopy.getNode(b1, taxaOrdering);
							MyNode v2 = tCopy.getNode(b2, taxaOrdering);
							if (!v1.equals(tCopy.getRoot())) {
								if (tCopy.pruneAndRegraft(v1, v2, taxaOrdering)) {
									cmpDistortionScoreRec(tCopy, taxaOrdering, A, B, new Integer(score + 1), k);
								}
							}
						}
					}
				}
			}
		}

	}

	private boolean containmentCheck(BitSet b1, BitSet b2) {
		BitSet b = (BitSet) b2.clone();
		b.and(b1);
		return b.equals(b2);
	}

	private boolean checkTree(MyTree t, ArrayList<String> taxaOrdering, BitSet A, BitSet B) {
		for (BitSet c : t.getAllCluster(taxaOrdering)) {
			if (c.equals(A) || c.equals(B))
				return true;
		}
		return false;
	}

}
