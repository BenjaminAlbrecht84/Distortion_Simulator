package util;

import java.util.ArrayList;

public class RandomBinaryTree_Generator {

	public static Object[] run(int numberOfTaxa) {
		ArrayList<MyNode> nodes = new ArrayList<MyNode>();
		ArrayList<String> taxaOrdering = new ArrayList<String>();
		int idCounter = 0;
		for (int i = 0; i < numberOfTaxa; i++) {
			String taxon = String.valueOf(idCounter++);
			nodes.add(new MyNode(taxon));
			taxaOrdering.add(taxon);
		}
		createInnerNodes(nodes, idCounter);
		Object[] res = { new MyTree(nodes.get(0)), taxaOrdering };
		return res;
	}

	private static void createInnerNodes(ArrayList<MyNode> nodes, int idCounter) {
		while (nodes.size() != 1) {
			MyNode v1 = getRandomNode(nodes);
			nodes.remove(v1);
			MyNode v2 = getRandomNode(nodes);
			nodes.remove(v2);

			MyNode v = new MyNode("");
			v.addChildren(v1);
			v.addChildren(v2);
			nodes.add(v);
		}
	}

	private static MyNode getRandomNode(ArrayList<MyNode> nodes) {
		double d = Math.random() * (nodes.size() - 1);
		int index = (int) Math.round(d);
		return nodes.get(index);
	}

}
