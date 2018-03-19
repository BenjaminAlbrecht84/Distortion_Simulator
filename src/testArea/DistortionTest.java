package testArea;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.Random;

import method.Distortion_Scorer;
import method.MyDistortion_Scorer;
import util.MyNode;
import util.MyTree;
import util.NaiveDistortion_Scorer;
import util.MyNewickParser;
import util.RandomBinaryTree_Generator;

public class DistortionTest {

	public static void main(String[] args) throws IOException {

		myTreeTest();
		// randomTest();
		// paperTest();
		// myTest();

	}

	private static void myTreeTest() {

		String s0 = "((NC_010159,NZ_CP009935)NODE_0000001,(((NZ_CP010023,NC_005810)NODE_0000012,(((NC_004088,NZ_CP006792)NODE_0000015,((NZ_CP009704,(NZ_CP009991,NC_017265)NODE_0000018)NODE_0000017,NC_008149)NODE_0000016)NODE_0000014,((NZ_CP009906,NC_008150)NODE_0000020,((NC_017154,(NC_017160,(((NZ_CP009723,NZ_CP009785)NODE_0000028,(NZ_CP016273,((NZ_CP009973,NC_003143)NODE_0000031,NZ_CP009844)NODE_0000030)NODE_0000029)NODE_0000027,(NZ_CP009492,(NC_017168,NZ_CP009840)NODE_0000026)NODE_0000025)NODE_0000024)NODE_0000023)NODE_0000022,NC_014029)NODE_0000021)NODE_0000019)NODE_0000013)NODE_0000011,(((NC_009381,NZ_CP009715)NODE_0000005,NZ_CP006758)NODE_0000004,((NZ_CP006754,((NZ_CP006751,NZ_CP006748)NODE_0000010,(NZ_CP006762,NZ_CP006783)NODE_0000009)NODE_0000008)NODE_0000007,NZ_CP010247)NODE_0000006)NODE_0000003)NODE_0000002):1.0;";
		String s1 = "((NZ_CP016273,((((NZ_CP009704,(NZ_CP009991,NC_017265)),NC_008149),(NC_004088,NZ_CP006792)),(NZ_CP009906,NC_008150),(NC_005810,NZ_CP010023),NC_014029,NC_017154)null,(((NC_009381,NZ_CP009715),NZ_CP006758),((NZ_CP006754,((NZ_CP006748,NZ_CP006751),(NZ_CP006762,NZ_CP006783))),NZ_CP010247)),((NC_017168,NZ_CP009840),NZ_CP009492),(NZ_CP009785,NZ_CP009723),(NC_010159,NZ_CP009935),NZ_CP009844,NC_017160,NC_003143)null,NZ_CP009973)null:0.0;";

		MyTree t0 = new MyNewickParser().run(s0);
		MyTree t1 = new MyNewickParser().run(s1);

		MyTree[] a = { t1 };
		ArrayList<MyTree> trees = new ArrayList<MyTree>(Arrays.asList(a));
		ArrayList<MyTree> allTrees = (ArrayList<MyTree>) trees.clone();
		ArrayList<String> taxa = new ArrayList<String>();
		allTrees.add(t0);
		for (MyTree t : allTrees) {
			for (MyNode v : t.getNodes()) {
				if (v.isLeaf()) {
					String label = v.getId();
					v.setID(label.split("\\|")[0]);
					if (!taxa.contains(label))
						taxa.add(label);
				} else
					v.setID("");
			}
		}
		Collections.sort(taxa);

		ArrayList<String> taxaOrdering = new ArrayList<String>();
		for (MyTree t : allTrees) {
			for (MyNode l : t.getLeafNodes()) {
				l.setID(String.valueOf(taxa.indexOf(l.getId())));
				String label = l.getId();
				if (!taxaOrdering.contains(label))
					taxaOrdering.add(label);
			}
		}
		// Collections.sort(taxaOrdering);
		Collections.sort(taxaOrdering, new TaxaComparator());

		System.out.println(t0.toNewickString());
		System.out.println(t1.toNewickString());

		MyDistortion_Scorer s = new MyDistortion_Scorer(t0, taxaOrdering);
		for (MyNode v : t1.getNodes()) {
			BitSet c = t1.getCluster(v, taxaOrdering);
			try {
				int score = s.run(c);
				System.out.println("> "+score+"\n"+c);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

	}

	private static void myTest() {

		String s = "(14:0.0,17:0.0,((9:0.0,(24:0.0,8:0.0):0.0):0.0,(10:0.0,(((((6:0.0,3:0.0):0.0,7:0.0):0.0,(19:0.0,(1:0.0,20:0.0):0.0):0.0):0.0,(((12:0.0,11:0.0):0.0,(15:0.0,(5:0.0,(2:0.0,0:0.0):0.0):0.0):0.0):0.0,(22:0.0,((13:0.0,16:0.0):0.0,(23:0.0,18:0.0):0.0):0.0):0.0):0.0):0.0,(4:0.0,21:0.0):0.0):0.0):0.0):0.0)null:0.0;";
		MyTree t = new MyNewickParser().run(s);

		ArrayList<Integer> taxa = new ArrayList<Integer>();
		for (MyNode l : t.getLeafNodes())
			taxa.add(Integer.parseInt(l.getId()));
		Collections.sort(taxa);
		ArrayList<String> taxaOrdering = new ArrayList<String>();
		for (int i : taxa)
			taxaOrdering.add(String.valueOf(i));

		String[] aSet = { "1", "3", "6", "7", "10", "19", "20" };
		BitSet A = new BitSet(taxaOrdering.size());
		for (String c : aSet)
			A.set(taxaOrdering.indexOf(c));

		BitSet B = new BitSet(taxaOrdering.size());
		B.set(0, taxaOrdering.size());
		B.xor(A);

		int counter = taxaOrdering.size();
		for (MyNode v : t.getNodes()) {
			if (!v.isLeaf())
				v.setID(String.valueOf(counter++));
		}

		System.out.println(t.toNewickString());
		int distortion = new Distortion_Scorer().run(t, A, B, taxaOrdering);
		System.out.println(distortion);

	}

	private static void randomTest() throws IOException {

		int[] rSPRMoves = { 20 };
		int[] taxaNumbers = { 10 };
		int iterations = 100;

		Random rand = new Random();

		for (int k : rSPRMoves) {

			for (int numberOfTaxa : taxaNumbers) {

				for (int i = 0; i < iterations; i++) {

					// creating random binary tree
					Object[] res = RandomBinaryTree_Generator.run(numberOfTaxa);
					MyTree t = (MyTree) res[0];
					ArrayList<String> taxaOrdering = (ArrayList<String>) res[1];

					// choosing split
					MyNode v = null;
					int lowerBound = (int) Math.round(new Double(numberOfTaxa) / 3.);
					int upperBound = numberOfTaxa - lowerBound;
					while (v == null || t.getCluster(v, taxaOrdering).cardinality() < lowerBound
							|| t.getCluster(v, taxaOrdering).cardinality() > upperBound) {
						ArrayList<MyNode> nodes = t.getNodes();
						v = nodes.get(rand.nextInt(nodes.size()));
					}

					int counter = 0;
					BitSet A = new BitSet(taxaOrdering.size());
					MyTree tCopy = new MyNewickParser().run(t.toNewickString());
					while (counter < k) {
						ArrayList<MyNode> nodes = tCopy.getNodes();
						MyNode v1 = nodes.get(rand.nextInt(nodes.size()));
						MyNode v2 = nodes.get(rand.nextInt(nodes.size()));
						BitSet b1 = tCopy.getCluster(v1, taxaOrdering);
						BitSet b2 = tCopy.getCluster(v2, taxaOrdering);
						if (!containmentCheck(b1, b2)) {
							tCopy.pruneAndRegraft(v1, v2, taxaOrdering);
							counter++;
						}
					}

					A = tCopy.getCluster(v, taxaOrdering);
					BitSet B = new BitSet(taxaOrdering.size());
					B.set(0, taxaOrdering.size());
					B.xor(A);

					int d1 = new NaiveDistortion_Scorer().run(tCopy, taxaOrdering, A, B);
					int d2 = new MyDistortion_Scorer(tCopy, taxaOrdering).run(A, B);

					if (d1 != d2)
						System.err.println("ERROR:");
					// if (d1 > 1 || d1 != d2) {
					System.out.println(t.toNewickString());
					System.out.println(tCopy.toNewickString());
					System.out.println(A + " " + B);
					System.out.println("->" + d1 + " " + d2);
					// }

				}

			}

		}

	}

	private static boolean containmentCheck(BitSet b1, BitSet b2) {
		BitSet b = (BitSet) b2.clone();
		b.and(b1);
		return b.equals(b2);
	}

	private static void cmpLeafSetRec(ArrayList<MyNode> leafSet, MyNode v, ArrayList<String> taxaOrdering, BitSet b) {
		if (v.isLeaf()) {
			leafSet.add(v);
			b.set(taxaOrdering.indexOf(v.getId()));
		}
		for (MyNode w : v.getChildren())
			cmpLeafSetRec(leafSet, w, taxaOrdering, b);
	}

	private static void paperTest() {

		String s = "((3,2,0,1),(6,7,8,(4,5)));";
		MyTree t = new MyNewickParser().run(s);

		String[] taxa = { "0", "1", "2", "3", "4", "5", "6", "7", "8" };
		ArrayList<String> taxaOrdering = new ArrayList<String>(Arrays.asList(taxa));

		String[] aSet = { "0", "1", "8", "5", "4" };
		BitSet A = new BitSet(taxaOrdering.size());
		for (String c : aSet)
			A.set(taxaOrdering.indexOf(c));

		String[] bSet = { "2", "3", "6", "7" };
		BitSet B = new BitSet(taxaOrdering.size());
		for (String c : bSet)
			B.set(taxaOrdering.indexOf(c));

		int counter = 0;
		for (MyNode v : t.getNodes()) {
			if (!v.isLeaf())
				v.setID(String.valueOf(counter++));
		}

		System.out.println(t.toNewickString());
		// int d1 = new Distortion_Scorer().run(t, A, B, taxaOrdering);
		int d2 = new NaiveDistortion_Scorer().run(t, taxaOrdering, A, B);
		// System.out.println(d1);

	}

	private static MyNode getRandomNode(ArrayList<MyNode> nodes) {
		double d = Math.random() * (nodes.size() - 1);
		int index = (int) Math.round(d);
		return nodes.get(index);
	}

}
