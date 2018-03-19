package simulation;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;

import util.MyNode;
import util.MyTree;
import util.MyNewickParser;

public class RAxML_TreeGenerator {

	private final static String RAxML_EXE = "/Users/Benjamin/Documents/workspace/Distortion_Simulator/tools/raxml";

	public ArrayList<MyTree> run(ArrayList<HashMap<MyNode, String>> setOfEvolvedSeqs) throws IOException {

		ArrayList<MyTree> trees = new ArrayList<MyTree>();
		for (HashMap<MyNode, String> MSA : setOfEvolvedSeqs) {

			File inputFile = File.createTempFile("raxmlInput", ".txt");

			// generating input file
			StringBuilder buf = new StringBuilder();
			buf.append(MSA.keySet().size() + "\n");
			buf.append(MSA.values().iterator().next().length() + "\n");
			for (MyNode v : MSA.keySet())
				buf.append(v.getId() + " " + MSA.get(v) + "\n");
			FileWriter fW = new FileWriter(inputFile);
			fW.write(buf.toString());
			fW.close();

			// running RAxML
			String cmd = RAxML_EXE + " -p 100 -T 8 -s " + inputFile.getAbsolutePath() + " -n txt" + " -m GTRCAT";
			executingCommand(cmd);

			// reading output file
			MyTree t = parsingOutputFile(new File("RAxML_result.txt"));
			trees.add(t);

			// deleting files
			inputFile.delete();
			String[] suffices = { "bestTree", "info", "log", "parsimonyTree", "result" };
			for (String s : suffices)
				new File("RAxML_" + s + ".txt").delete();

		}

		return trees;

	}

	private MyTree parsingOutputFile(File f) throws IOException {

		BufferedReader buf = new BufferedReader(new FileReader(f));
		String l = buf.readLine();
		MyTree t = new MyNewickParser().run(l);
		buf.close();

		return t;

	}

	private int executingCommand(String command) {
		try {

			// System.out.println("OUTPUT>Executing " + command);
			Runtime rt = Runtime.getRuntime();
			Process proc = rt.exec(command);

			// checking error messages
			StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "ERROR");

			// checking error messages
			StreamGobbler outputGobbler = new StreamGobbler(proc.getInputStream(), "OUTPUT");

			errorGobbler.start();
			outputGobbler.start();
			int exitVal = proc.waitFor();

			return exitVal;

		} catch (Exception e) {
			e.printStackTrace();
		}

		return 1;
	}

	class StreamGobbler extends Thread {
		InputStream is;
		String type;

		StreamGobbler(InputStream is, String type) {
			this.is = is;
			this.type = type;
		}

		public void run() {
			try {
				InputStreamReader isr = new InputStreamReader(is);
				BufferedReader br = new BufferedReader(isr);
				String line = null;
				while ((line = br.readLine()) != null) {
					if (type.equals("ERROR"))
						System.out.println(type + ">" + line);
					// if (type.equals("OUTPUT"))
					// System.out.println(type + ">" + line);
				}
			} catch (IOException ioe) {
				ioe.printStackTrace();
			}
		}
	}

}
