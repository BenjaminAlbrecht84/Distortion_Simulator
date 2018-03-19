package util;

import java.util.ArrayList;

public class MyNode {

	private String id = "";
	private double branchLength = 1.;
	private ArrayList<MyNode> children = new ArrayList<MyNode>();
	private MyNode parent;

	private Object[] distortionOperation;

	public MyNode() {
	}

	public MyNode(String id) {
		this.id = id;
	}

	public void addChildren(MyNode v) {
		if (!children.contains(v)) {
			children.add(v);
			v.setParent(this);
		}
	}

	public String toNewick(String newickString) {
		if (children.isEmpty() && parent == null)
			return "(" + newickLabel() + ")";
		else if (children.isEmpty()) {
			return newickLabel();
		} else {
			String subString = "(";
			for (MyNode c : children) {
				if (children.indexOf(c) == children.size() - 1)
					subString = subString.concat(c.toNewick(newickString) + ")");
				else
					subString = subString.concat(c.toNewick(newickString) + ",");
			}
			return subString + newickLabel();
		}
	}

	public void removeChild(MyNode c) {
		children.remove(c);
		c.setParent(null);
	}

	public MyNode copy() {
		MyNode vCopy = new MyNode();
		vCopy.setID(this.id);
		vCopy.setBranchLength(this.getBranchLength());
		return vCopy;
	}

	private String newickLabel() {
		return id + ":" + branchLength;
	}

	public ArrayList<MyNode> getChildren() {
		return children;
	}

	public void setChildren(ArrayList<MyNode> children) {
		this.children = children;
	}

	public MyNode getParent() {
		return parent;
	}

	public void setParent(MyNode parent) {
		this.parent = parent;
	}

	public String getId() {
		return id;
	}

	public boolean isLeaf() {
		return children.isEmpty();
	}

	public double getBranchLength() {
		return branchLength;
	}

	public void setBranchLength(double branchLength) {
		this.branchLength = branchLength;
	}

	public void setID(String id) {
		this.id = id;
	}

	public Object[] getDistortionOperation() {
		return distortionOperation;
	}

	public void setDistortionOperation(Object[] distortionOperation) {
		this.distortionOperation = distortionOperation;
	}

}
