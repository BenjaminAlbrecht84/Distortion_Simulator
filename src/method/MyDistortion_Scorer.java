package method;

import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;

import util.MyNode;
import util.MyTree;

public class MyDistortion_Scorer {

	private int memoryCounter = 0;
	private HashMap<BitSet, Integer> distortionMemory = new HashMap<BitSet, Integer>();
	private MyTree t;
	private ArrayList<String> taxaOrdering;

	public MyDistortion_Scorer(MyTree t, ArrayList<String> taxaOrdering) {
		this.t = t;
		this.taxaOrdering = taxaOrdering;
	}

	public int run(BitSet A) throws IOException {

		BitSet B = new BitSet(taxaOrdering.size());
		B.set(0, taxaOrdering.size());
		B.xor(A);

		if (distortionMemory.containsKey(A)) {
			memoryCounter++;
			return distortionMemory.get(A);
		}
		if (distortionMemory.containsKey(B)) {
			memoryCounter++;
			return distortionMemory.get(B);
		}

		return run(A, B);
	}

	public int run(BitSet A, BitSet B) throws IOException {

		HashMap<MyNode, Integer> scoreA = new HashMap<MyNode, Integer>();
		HashMap<MyNode, Integer> scoreB = new HashMap<MyNode, Integer>();
		ArrayList<MyNode> leaves = t.getLeafNodes();

		// checking overlap of splits with leaf-set
		BitSet leafSet = new BitSet(taxaOrdering.size());
		for (MyNode l : leaves)
			leafSet.set(taxaOrdering.indexOf(l.getId()));
		BitSet b = (BitSet) leafSet.clone();
		b.and(A);
		if (b.cardinality() <= 1)
			return 0;
		b = (BitSet) leafSet.clone();
		b.and(B);
		if (b.cardinality() <= 1)
			return 0;

		// initializing A- & B-scores for each leaf
		for (MyNode leaf : leaves) {

			// looking into leaf-sets
			int index = taxaOrdering.indexOf(leaf.getId());
			boolean hasA = A.get(index), hasB = B.get(index);

			// setting A- & B-values
			int aValue = 0, bValue = 0;
			if (hasA && !hasB)
				bValue = 1;
			else if (!hasA && hasB)
				aValue = 1;
			else { // there is something wrong with the input
				if (hasA && hasB)
					throw new IOException("ERROR: Taxon t=" + leaf.getId() + " is present in both splits.");
				if (!hasA && !hasB)
					throw new IOException("ERROR: Taxon t=" + leaf.getId() + " not present in both splits.");
				return -1;
			}

			scoreA.put(leaf, aValue);
			scoreB.put(leaf, bValue);

		}

		// computing distortion score
		computeScoresRec(t.getRoot(), scoreA, scoreB);
		int rootScoreA = scoreA.get(t.getRoot());
		int rootScoreB = scoreB.get(t.getRoot());

		int distortion = Math.min(rootScoreA, rootScoreB) - 1;

		distortionMemory.put(A, distortion);
		distortionMemory.put(B, distortion);
		return distortion;

	}

	private void computeScoresRec(MyNode v, HashMap<MyNode, Integer> scoreA, HashMap<MyNode, Integer> scoreB) {

		boolean isABeatenByB = false, isBBeatenByA = false;
		int countA = 0, countB = 0;

		if (!v.isLeaf()) {

			for (MyNode w : v.getChildren()) {
				computeScoresRec(w, scoreA, scoreB);
				if (scoreA.get(w) < scoreB.get(w)) {
					isBBeatenByA = true;
					countB += scoreA.get(w);
				} else
					countB += scoreB.get(w);
				if (scoreB.get(w) < scoreA.get(w)) {
					isABeatenByB = true;
					countA += scoreB.get(w);
				} else
					countA += scoreA.get(w);
			}

			if (isBBeatenByA)
				countB++;
			if (isABeatenByB)
				countA++;

			scoreB.put(v, countB);
			scoreA.put(v, countA);
			// System.out.println(v.getId() + " :" + scoreA.get(v) + "/" + scoreB.get(v));

		}

	}

	public int getMemoryCounter() {
		return memoryCounter;
	}

}
