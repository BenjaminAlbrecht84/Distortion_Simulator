package util;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.Random;

public class MyTree {

	private MyNode root, lastRoot;
	private HashMap<BitSet, MyNode> cluster2node = new HashMap<BitSet, MyNode>();

	public MyTree(MyNode root) {
		this.root = root;
	}

	public void attachNodeAbove(MyNode v, MyNode w, double[] branchLengths) {
		MyNode x = new MyNode();
		MyNode p = w.getParent();
		if (p != null) {
			p.removeChild(w);
			p.addChildren(x);
		}
		x.addChildren(w);
		x.addChildren(v);
		if (p == null)
			root = x;
		if (branchLengths != null) {
			x.setBranchLength(branchLengths[0]);
			w.setBranchLength(branchLengths[1]);
		}
	}

	public MyTree extractSubtree(MyNode v) {
		MyNode vCopy = v.copy();
		MyTree tCopy = new MyTree(vCopy);
		extractSubtreeRec(vCopy, v);
		return tCopy;
	}

	private void extractSubtreeRec(MyNode vCopy, MyNode v) {
		for (MyNode w : v.getChildren()) {
			MyNode wCopy = w.copy();
			vCopy.addChildren(wCopy);
			extractSubtreeRec(wCopy, w);
		}
	}

	public void assingBranchLength(int maxLength, int numOfLeaves) {
		double h = Math.sqrt(new Double(numOfLeaves));
		assingBranchLengthRec(root, maxLength, h);
	}

	public ArrayList<BitSet> getAllCluster(ArrayList<String> taxaOrdering) {
		cluster2node = new HashMap<BitSet, MyNode>();
		cmpAllClustersRec(root, taxaOrdering);
		return new ArrayList<BitSet>(cluster2node.keySet());
	}

	public MyNode getNode(BitSet b, ArrayList<String> taxaOrdering) {
		cluster2node = new HashMap<BitSet, MyNode>();
		cmpAllClustersRec(root, taxaOrdering);
		return cluster2node.get(b);
	}

	public BitSet getCluster(MyNode v, ArrayList<String> taxaOrdering) {
		cluster2node = new HashMap<BitSet, MyNode>();
		return cmpAllClustersRec(v, taxaOrdering);
	}

	public MyNode getLCA(BitSet cluster, ArrayList<String> taxaOrdering) {
		cluster2node = new HashMap<BitSet, MyNode>();
		cmpAllClustersRec(root, taxaOrdering);
		BitSet lcaCluster = null;
		for (BitSet c : cluster2node.keySet()) {
			BitSet b = (BitSet) c.clone();
			b.or(cluster);
			if (b.equals(c) && (lcaCluster == null || c.cardinality() < lcaCluster.cardinality()))
				lcaCluster = c;
		}
		return cluster2node.get(lcaCluster);
	}

	public double getAvgEdgeLength() {
		ArrayList<MyNode> nodes = getNodes();
		double sum = 0;
		for (MyNode v : nodes)
			sum += v.getBranchLength();
		return sum / (new Double(nodes.size() - 1.));
	}

	public ArrayList<MyNode> getLeafNodes() {
		ArrayList<MyNode> leaves = new ArrayList<MyNode>();
		getLevesRec(root, leaves);
		return leaves;
	}

	private void getLevesRec(MyNode v, ArrayList<MyNode> leaves) {
		if (v.isLeaf())
			leaves.add(v);
		for (MyNode w : v.getChildren())
			getLevesRec(w, leaves);
	}

	public ArrayList<MyNode> getNodes() {
		ArrayList<MyNode> nodes = new ArrayList<MyNode>();
		getNodesRec(root, nodes);
		return nodes;
	}

	private void getNodesRec(MyNode v, ArrayList<MyNode> nodes) {
		nodes.add(v);
		for (MyNode c : v.getChildren())
			getNodesRec(c, nodes);
	}

	private BitSet cmpAllClustersRec(MyNode v, ArrayList<String> taxaOrdering) {
		BitSet b = new BitSet(taxaOrdering.size());
		if (v.isLeaf()) {
			if (!v.getId().isEmpty() && taxaOrdering.contains(v.getId()))
				b.set(taxaOrdering.indexOf(v.getId()));
		} else {
			for (MyNode c : v.getChildren())
				b.or(cmpAllClustersRec(c, taxaOrdering));
		}
		if (!isIsolatedNode(v) && !cluster2node.containsKey(b))
			cluster2node.put(b, v);
		return b;
	}

	private boolean isIsolatedNode(MyNode v) {
		if (v.getChildren().size() == 1 && v.getParent() != null)
			return true;
		return false;
	}

	public boolean pruneAndRegraft(MyNode v1, MyNode v2, ArrayList<String> taxaOrdering) {

		if (v1.getParent() != null && v1.getParent().equals(v2) && v2.getChildren().size() == 2)
			return false;

		MyNode p1 = v1.getParent();
		p1.removeChild(v1);

		supressNode(p1, taxaOrdering);

		if (v2.getParent() != null) {
			MyNode p2 = v2.getParent();
			p2.removeChild(v2);
			MyNode x = new MyNode();
			p2.addChildren(x);
			x.addChildren(v1);
			x.addChildren(v2);
		} else {
			root = new MyNode();
			root.addChildren(v2);
			root.addChildren(v1);
		}

		return true;

	}

	public void removeNode(MyNode v, boolean supressNode, ArrayList<String> taxaOrdering) {
		if (v.getParent() != null) {
			MyNode p = v.getParent();
			p.removeChild(v);
			if (supressNode)
				supressNode(p, taxaOrdering);
		} else
			root = null;
	}

	public boolean supressNode(MyNode v, ArrayList<String> taxaOrdering) {

		if (v.getChildren().size() == 1 && v.getParent() != null) {

			MyNode p = v.getParent();
			MyNode c = v.getChildren().get(0);
			v.removeChild(c);
			p.removeChild(v);
			p.addChildren(c);
			c.setBranchLength(c.getBranchLength() + v.getBranchLength());

			return true;

		} else if (v.getChildren().size() == 1) {

			MyNode c = v.getChildren().get(0);
			c.setParent(null);
			root = c;

			return true;

		}
		return false;
	}

	private void assingBranchLengthRec(MyNode v, int totalLength, double h) {

		int rootDist = 0;
		MyNode p = v.getParent();
		while (p != null) {
			rootDist += p.getBranchLength();
			p = p.getParent();
		}

		if (v.getChildren().isEmpty())
			v.setBranchLength(totalLength - rootDist);
		else {
			if (v.getParent() != null) {
				Random rand = new Random();
				int randBranchLength = rand.nextInt(totalLength - rootDist);
				int max = (int) Math.round((1. / h) * totalLength);
				int branchLength = randBranchLength > max ? max : randBranchLength;
				v.setBranchLength(branchLength);
			}
			for (MyNode w : v.getChildren())
				assingBranchLengthRec(w, totalLength, h);
		}

	}

	public MyTree copy() {

		// copying tree
		HashMap<MyNode, MyNode> node2copy = new HashMap<MyNode, MyNode>();
		MyNode rootCopy = new MyNode();
		MyTree tCopy = new MyTree(rootCopy);
		copyRec(root, rootCopy, node2copy);

		// copying distortion links
		for (MyNode v : getNodes()) {
			MyNode vCopy = node2copy.get(v);
			if (v.getDistortionOperation() != null) {
				MyNode link = node2copy.get(v.getDistortionOperation()[0]);
				int order = (int) v.getDistortionOperation()[1];
				Double distance = v.getDistortionOperation()[2] == null ? null : (Double) v.getDistortionOperation()[2];
				Object[] distortionOperation = { link, order, distance };
				vCopy.setDistortionOperation(distortionOperation);
			}
		}
		if (lastRoot != null)
			tCopy.setLastRoot(node2copy.get(lastRoot));

		return tCopy;
	}

	private void copyRec(MyNode v, MyNode vCopy, HashMap<MyNode, MyNode> node2copy) {
		vCopy.setID(v.getId());
		vCopy.setBranchLength(v.getBranchLength());
		node2copy.put(v, vCopy);
		for (MyNode w : v.getChildren()) {
			MyNode wCopy = new MyNode();
			vCopy.addChildren(wCopy);
			copyRec(w, wCopy, node2copy);
		}
	}

	public MyNode getRoot() {
		return root;
	}

	public String toNewickString() {
		return root.toNewick("") + ";";
	}

	public void setRoot(MyNode r) {
		this.root = r;
	}

	public MyNode getLastRoot() {
		return lastRoot;
	}

	public void setLastRoot(MyNode lastRoot) {
		this.lastRoot = lastRoot;
	}

	public ArrayList<String> getTaxa() {
		ArrayList<String> taxa = new ArrayList<String>();
		for (MyNode l : this.getLeafNodes()) {
			if (!l.getId().isEmpty())
				taxa.add(l.getId());
		}
		return taxa;
	}

}
