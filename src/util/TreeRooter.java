package util;

import java.util.ArrayList;
import java.util.BitSet;

public class TreeRooter {

	public static void run(MyTree t, MyNode newRoot, ArrayList<String> taxaOrdering) {

		if (newRoot.equals(t.getRoot()))
			return;

		ArrayList<MyNode> rootPath = new ArrayList<MyNode>();
		MyNode p = newRoot.getParent();
		while (p != null) {
			rootPath.add(p);
			p = p.getParent();
		}

		// re-rooting tree
		for (int i = 0; i < rootPath.size() - 1; i++) {
			rootPath.get(i + 1).removeChild(rootPath.get(i));
			rootPath.get(i).addChildren(rootPath.get(i + 1));
			if (i > 0)
				rootPath.get(i).setParent(rootPath.get(i - 1));
		}
		newRoot.getParent().removeChild(newRoot);

		// building left and right tree
		MyTree t2 = new MyTree(rootPath.get(0));
		MyNode lastRoot = t.getRoot().getChildren().size() == 1 ? t.getRoot().getChildren().get(0) : t.getRoot();
		t2.supressNode(t.getRoot(), taxaOrdering);
		MyTree t3 = new MyTree(newRoot);

		// generating re-rooted tree
		MyNode root = new MyNode();
		MyTree t4 = new MyTree(root);
		root.addChildren(t2.getRoot());
		root.addChildren(t3.getRoot());
		t.setLastRoot(lastRoot);
		t.setRoot(t4.getRoot());

	}

	public static void run(MyTree t, BitSet cluster, ArrayList<String> taxaOrdering) {
		rerootTree(t, cluster, taxaOrdering);
	}

	public static void run(MyTree t2, MyTree t1, ArrayList<String> taxaOrdering) {
		for (MyNode w : t1.getRoot().getChildren()) {
			BitSet cluster = t1.getCluster(w, taxaOrdering);
			if (t2.getNode(cluster, taxaOrdering) == null) {
				rerootTree(t2, cluster, taxaOrdering);
				break;
			}
		}
	}

	private static void rerootTree(MyTree t, BitSet cluster, ArrayList<String> taxaOrdering) {

		// computing root cluster
		BitSet rootSet = new BitSet(taxaOrdering.size());
		rootSet.set(0, taxaOrdering.size());
		rootSet.xor(cluster);

		// collecting information of original tree
		MyNode lca = t.getLCA(rootSet, taxaOrdering);
		if (lca.equals(t.getRoot()))
			return;
		resolveTree(t, rootSet, taxaOrdering);
		ArrayList<MyNode> rootPath = new ArrayList<MyNode>();

		MyNode p = lca.getParent();
		while (p != null) {
			rootPath.add(p);
			p = p.getParent();
		}

		// re-rooting tree
		for (int i = 0; i < rootPath.size() - 1; i++) {
			rootPath.get(i + 1).removeChild(rootPath.get(i));
			rootPath.get(i).addChildren(rootPath.get(i + 1));
			if (i > 0)
				rootPath.get(i).setParent(rootPath.get(i - 1));
		}
		lca.getParent().removeChild(lca);

		// building left and right tree
		MyTree t2 = new MyTree(rootPath.get(0));
		MyNode lastRoot = t.getRoot().getChildren().size() == 1 ? t.getRoot().getChildren().get(0) : t.getRoot();
		t2.supressNode(t.getRoot(), taxaOrdering);
		MyTree t3 = new MyTree(lca);

		// generating re-rooted tree
		MyNode root = new MyNode();
		MyTree t4 = new MyTree(root);
		root.addChildren(t2.getRoot());
		root.addChildren(t3.getRoot());
		t.setLastRoot(lastRoot);
		t.setRoot(t4.getRoot());

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
			rerootTree(t, cluster, taxaOrdering);

	}

}
