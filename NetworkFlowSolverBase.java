package testing;

/**
 * @author William Fiset, william.alexandre.fiset@gmail.com
 **/
//package com.williamfiset.algorithms.graphtheory.networkflow;

import static java.lang.Math.min;

import java.util.ArrayList;
import java.util.List;

public abstract class NetworkFlowSolverBase {

	// To avoid overflow, set infinity to a value less than Long.MAX_VALUE;
	protected static final long INF = Long.MAX_VALUE / 2;

	public static class Edge {
		public int from, to;
		public Edge residual;
		public long flow, cost;
		public final long capacity, originalCost;

		public Edge(int from, int to, long capacity) {
			this(from, to, capacity, 0 /* unused */);
		}

		public Edge(int from, int to, long capacity, long cost) {
			this.from = from;
			this.to = to;
			this.capacity = capacity;
			this.originalCost = this.cost = cost;
		}

		public boolean isResidual() {
			return capacity == 0;
		}

		public long remainingCapacity() {
			return capacity - flow;
		}

		public void augment(long bottleNeck) {
			flow += bottleNeck;
			residual.flow -= bottleNeck;
		}

		public String toString(int s, int t) {
//			String u = (from == s) ? "s" : ((from == t) ? "t" : String.valueOf(from));
//			String v = (to == s) ? "s" : ((to == t) ? "t" : String.valueOf(to));
			String u = String.valueOf(from);
			String v =  String.valueOf(to);
			return String.format("Edge %s -> %s | flow = %d | capacity = %d | is residual: %s", 
					u, v, flow, capacity, isResidual());
			
		}
		
		public boolean isMinCut() {
			return capacity == flow;
		}
	}

	// Inputs: n = number of nodes, s = source, t = sink
	protected final int n, s, t;
	protected final int bigN, bigS, bigT;

	protected long maxFlow;
	protected long bigMaxFlow;
	protected long minCost; //not relevant

	//nodes left of cut
	protected boolean[] leftOfCut;
	protected boolean[] bigLeftOfCut;

	//array of adjacency lists, representing the graph
	protected List<Edge>[] graph;
	protected List<Edge>[] bigGraph;

	protected ArrayList<Integer> originalCutNodes;
	protected ArrayList<Edge> cutEdges;

	// 'visited' and 'visitedToken' are variables used for graph sub-routines to 
	// track whether a node has been visited or not. In particular, node 'i' was 
	// recently visited if visited[i] == visitedToken is true. This is handy 
	// because to mark all nodes as unvisited simply increment the visitedToken.
	private int visitedToken = 1;
	private int bigVisitedToken = 1;
	private int[] visited;
	private int[] bigVisited;


	// Indicates whether the network flow algorithm has ran. We should not need to
	// run the solver multiple times, because it always yields the same result.
	private boolean solvedBFS;	
	private boolean bigSolvedBFS;


	/**
	 * Creates an instance of a flow network solver. Use the {@link #addEdge}
	 * method to add edges to the graph.
	 *
	 * @param n - The number of nodes in the graph including source and sink nodes.
	 * @param s - The index of the source node, 0 <= s < n
	 * @param t - The index of the sink node, 0 <= t < n, t != s
	 */
	public NetworkFlowSolverBase(int n, int s, int t) {
		this.n = n; this.s = s; this.t = t; 
		this.bigN = n * 2; this.bigS = s * 2; this.bigT = t * 2; 
		initializeGraph();
		leftOfCut = new boolean[n];
		visited = new int[n];
		bigLeftOfCut = new boolean[bigN];
		bigVisited = new int[bigN];
		originalCutNodes = new ArrayList<Integer>();
		cutEdges =  new ArrayList<Edge>();
	}

	// Construct an empty graph with n nodes including the source and sink nodes.
	private void initializeGraph() {
		graph = new List[n];
		for (int i = 0; i < graph.length; i++)
			graph[i] = new ArrayList<Edge>();
	}

	/**
	 * Adds a directed edge (and residual edge) to the flow graph.
	 *
	 * @param from     - The index of the node the directed edge starts at.
	 * @param to       - The index of the node the directed edge ends at.
	 * @param capacity - The capacity of the edge.
	 */
	public void addEdge(int from, int to, long capacity) {
		if (capacity < 0) throw new IllegalArgumentException("Capacity < 0");
		Edge e1 = new Edge(from, to, capacity);
		Edge e2 = new Edge(to, from, 0);
		e1.residual = e2;
		e2.residual = e1;
		graph[from].add(e1);
		graph[to].add(e2);
	}
	
	/*
	public void deepGraphBackup(List<Edge>[] bigGraph) {
		for (int i = 1; i < bigGraph.length; i++) 
			for (Edge e : bigGraph[i]) {
				Edge temp = new Edge(e.from, e.to, e.capacity);
				bigGraphBackup[i].add(temp);	
			}
				 					
	}
	
	public void deepGraphRestore() {
		for (int i = 1; i < bigGraphBackup.length; i++) 
			for (Edge e : bigGraphBackup[i]) {
				Edge temp = new Edge(e.from, e.to, e.capacity);
				bigGraph[i].add(temp);	
			}				 					
	}*/
	
	//create in solver a big graph based on conversion of real graph - each vertex becomes two vertices and edge.
	//create edge from each vertex, make it minimal capacity (1) so it will be selected as min cut
	//so we can know which nodes in original graph are the cut nodes
	public void initBigGraph() {
		bigGraph = new List[bigN];
		for (int i = 0; i < bigGraph.length; i++) {
			bigGraph[i] = new ArrayList<Edge>();
		}
		//add two vertices with edge per each original vertex
		for (int i = 0; i < graph.length; i++) {
			int newCapacity = 1;
			if(i == s || i == t)
				newCapacity = bigGraph.length * 2; //if its sink or source node we make its edge capacity large
			int newFrom = (i * 2);
			int newTo = (i * 2) + 1;
			Edge e1 = new Edge(newFrom, newTo, newCapacity); 
			Edge e2 = new Edge(newTo, newFrom, 0);
			e1.residual = e2;
			e2.residual = e1;
			bigGraph[newFrom].add(e1);
			bigGraph[newTo].add(e2);
			//convert existing edges
			for(Edge edge : graph[i]) {	
				if(!edge.isResidual()) {
					//make original edges capacity very big so NOT selected as min cut of bigGraph
					Edge updatedEdge = new Edge((edge.from * 2) + 1, edge.to * 2, bigGraph.length * 2); 
					Edge updatedEdgeResidual = new Edge(edge.to * 2, (edge.from * 2) + 1, 0);
					updatedEdge.residual = updatedEdgeResidual;
					updatedEdgeResidual.residual = updatedEdge;
					bigGraph[(edge.from * 2) + 1].add(updatedEdge);
					bigGraph[edge.to * 2].add(updatedEdgeResidual);
				}
			}			
		}
		
		
		//print big graph
//		for (int i = 0; i < bigGraph.length; i++) {
//			System.out.println("Big graph node i = " + i);
//			for (Edge e : bigGraph[i]) {
//				System.out.println(e.toString(e.from, e.to));				
//			}					
//		}
				
	}
	
	/** Cost variant of {@link #addEdge(int, int, int)} for min-cost max-flow */
	public void addEdge(int from, int to, long capacity, long cost) {
		Edge e1 = new Edge(from, to, capacity, cost);
		Edge e2 = new Edge(to, from, 0, -cost);
		e1.residual = e2;
		e2.residual = e1;
		graph[from].add(e1);
		graph[to].add(e2);
	}

	// Marks node 'i' as visited.
	public void visit(int i) {
		visited[i] = visitedToken;
	}
	
	public void bigVisit(int i) {
		bigVisited[i] = bigVisitedToken;
	}

	// Returns whether or not node 'i' has been visited.
	public boolean visited(int i) {
		return visited[i] == visitedToken;
	}
	
	public boolean bigVisited(int i) {
		return bigVisited[i] == bigVisitedToken;
	}

	// Resets all nodes as unvisited. This is especially useful to do
	// between iterations of finding augmenting paths, O(1)
	public void markAllNodesAsUnvisited() {
		visitedToken++;
	}
	
	public void bigMarkAllNodesAsUnvisited() {
		bigVisitedToken++;
	}

	/**
	 * Returns the graph after the solver has been executed. This allow you to
	 * inspect the {@link Edge#flow} compared to the {@link Edge#capacity} in 
	 * each edge. This is useful if you want to figure out which edges were 
	 * used during the max flow.
	 */
	public List<Edge>[] getGraph(boolean big) {
		execute(big);
		return big ? bigGraph : graph;
	}

	// Returns the maximum flow from the source to the sink.
	public long getMaxFlow(boolean big) {
		execute(big);
		return big ? bigMaxFlow : maxFlow;
	}

	// Returns the min cost from the source to the sink.
	// NOTE: This method only applies to min-cost max-flow algorithms.
	public long getMinCost(boolean big) {
		execute(big);
		return minCost;
	}

	// Returns the min-cut of this flow network in which the nodes on the "left side"
	// of the cut with the source are marked as true and those on the "right side" 
	// of the cut with the sink are marked as false.
	public boolean[] getMinCut(boolean big) {
		execute(big);
		return big ? bigLeftOfCut : leftOfCut;
	}

	// Wrapper method that ensures we only call solve() once
	private void execute(boolean big) {
		if (solvedBFS && !big || bigSolvedBFS && big) return;
		if (!big) {
			solvedBFS = true;
			solveBFS();	
			return;
		}
		//big
		initBigGraph(); //build big graph based on original graph
		bigSolvedBFS = true;
		bigSolveBFS();
	}


	// Method to implement which solves the network flow problem.
	public abstract void solveBFS();
	
	public abstract void bigSolveBFS();


}