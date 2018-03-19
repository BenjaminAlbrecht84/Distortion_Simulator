package util;

import java.util.ArrayList;
import java.util.BitSet;

public class Splits_Calculator {

	public static ArrayList<BitSet[]> run(MyTree t, ArrayList<String> taxaOrdering) {

		ArrayList<BitSet[]> splits = new ArrayList<BitSet[]>();
		computeSplitsRec(t.getRoot(), splits, taxaOrdering);
		return splits;

	}

	private static BitSet computeSplitsRec(MyNode v, ArrayList<BitSet[]> splits, ArrayList<String> taxaOrdering) {
		BitSet b = new BitSet(taxaOrdering.size());
		if (v.isLeaf())
			b.set(taxaOrdering.indexOf(v.getId()));
		else {
			for (MyNode w : v.getChildren())
				b.or(computeSplitsRec(w, splits, taxaOrdering));
		}
		BitSet[] split = createSplit(b, taxaOrdering);
		splits.add(split);
		return b;
	}

	private static BitSet[] createSplit(BitSet b1, ArrayList<String> taxaOrdering) {
		BitSet b2 = new BitSet(taxaOrdering.size());
		b2.set(0, taxaOrdering.size());
		b2.xor(b1);
		BitSet[] split = { b1, b2 };
		return split;
	}

}
