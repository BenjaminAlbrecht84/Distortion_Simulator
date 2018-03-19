package testArea;

import java.util.Comparator;

public class TaxaComparator implements Comparator<String> {

	@Override
	public int compare(String s1, String s2) {
		int i1 = Integer.parseInt(s1);
		int i2 = Integer.parseInt(s2);
		return Integer.compare(i1, i2);
	}

}
