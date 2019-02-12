package testing;
/**
 * An implementation of the Edmonds-Karp algorithm which is essentially
 * Ford-Fulkerson with a BFS as a method of finding augmenting paths. 
 * This Edmonds-Karp algorithm will allow you to find the max flow through
 * a directed graph and the min cut as a byproduct. 
 *
 * Time Complexity: O(VE^2)
 * 
 * @author William Fiset, william.alexandre.fiset@gmail.com
 **/
//package com.williamfiset.algorithms.graphtheory.networkflow;

import static java.lang.Math.min;

import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.LineNumberReader;
import java.io.StreamTokenizer;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.TimeUnit;

import testing.NetworkFlowSolverBase.Edge;

public class EdmondsKarpAdjacencyList extends NetworkFlowSolverBase {

	/**
	 * Creates an instance of a flow network solver. Use the {@link #addEdge(int, int, int)}
	 * method to add edges to the graph.
	 *
	 * @param n - The number of nodes in the graph including source and sink nodes.
	 * @param s - The index of the source node, 0 <= s < n
	 * @param t - The index of the sink node, 0 <= t < n, t != s
	 */
	public EdmondsKarpAdjacencyList(int n, int s, int t) {
		super(n, s, t);
	}

	// Run Edmonds-Karp and compute the max flow from the source to the sink node.
	@Override
	public void solveBFS() {
		long flow;
		do {
			markAllNodesAsUnvisited();
			flow = bfs(maxFlow);
			maxFlow += flow;
		} while (flow != 0);

		for(int i = 0; i < n; i++)
			if (visited(i))
				leftOfCut[i] = true;
	}
	
	@Override
	public void bigSolveBFS() {
		long flow;
		do {
			bigMarkAllNodesAsUnvisited();
			flow = bigBfs(bigMaxFlow);
			bigMaxFlow += flow;
		} while (flow != 0);

		for(int i = 0; i < bigN; i++) {
			if (bigVisited(i)) {
				bigLeftOfCut[i] = true;
				for(Edge edge : cutEdges) {
					if(edge.to == i)
						cutEdges.remove(edge);
				}
			}
		}
		
		for(Edge edge : cutEdges) {
			originalCutNodes.add(edge.from / 2);
		}
		
		
	}
	
	private long floww() {
		long flow = 0;
		for (Edge e : cutEdges)
			flow += e.capacity;
		return flow;	
	}

	private long bfs(long CurrentFlow) {
		Edge[] prev = new Edge[n];
		cutEdges.clear();
		originalCutNodes.clear();
		// The queue can be optimized to use a faster queue
		Queue<Integer> q = new ArrayDeque<>(n);
		visit(s);
		q.offer(s);

		// Perform BFS from source to sink
		while(!q.isEmpty()) {
			int node = q.poll();
			if (node == t) break;

			for (Edge edge : graph[node]) {
				long cap = edge.remainingCapacity();
				if (cap > 0 && !visited(edge.to)) {
					visit(edge.to);
					prev[edge.to] = edge;
					q.offer(edge.to);
				}
				if (cap == 0 && !edge.isResidual() && floww() < CurrentFlow) {
					cutEdges.add(edge);
					if(!originalCutNodes.contains(edge.from))
						originalCutNodes.add(edge.from);					
				}								
			}
		}

		// Sink not reachable!
		if (prev[t] == null) return 0;

		long bottleNeck = Long.MAX_VALUE;

		// Find augmented path and bottle neck
		for(Edge edge = prev[t]; edge != null; edge = prev[edge.from])
			bottleNeck = min(bottleNeck, edge.remainingCapacity());

		// Retrace augmented path and update flow values.
		for(Edge edge = prev[t]; edge != null; edge = prev[edge.from])
			edge.augment(bottleNeck);

		// Return bottleneck flow
		return bottleNeck;
	}

	private long bigBfs(long CurrentFlow) {
		Edge[] prev = new Edge[bigN];
		cutEdges.clear();
		// The queue can be optimized to use a faster queue
		Queue<Integer> q = new ArrayDeque<>(bigN);
		bigVisit(bigS);
		q.offer(bigS);
		
		System.out.println("START bfsssssssssss, current flow: " + CurrentFlow );
		// Perform BFS from source to sink
		while(!q.isEmpty()) {
			int node = q.poll();
			if (node == bigT) break;
			System.out.println("node " + node );
			for (Edge edge : bigGraph[node]) {
				long cap = edge.remainingCapacity();
				if (cap > 0 && !bigVisited(edge.to)) {
					bigVisit(edge.to);
					prev[edge.to] = edge;
					q.offer(edge.to);
				}
				System.out.println("edge " + edge.toString(edge.from, edge.to) );
				if (cap == 0 && !edge.isResidual() && floww() <= CurrentFlow) {
					Edge x = prev[edge.to];
					while(x != null) {
						if(cutEdges.contains(x))
							cutEdges.remove(x);
						x = prev[x.to];
					}
					cutEdges.add(edge);				
				}								
			}
		}
		
		// Sink not reachable!
		if (prev[bigT] == null) return 0;

		long bottleNeck = Long.MAX_VALUE;

		// Find augmented path and bottle neck
		for(Edge edge = prev[bigT]; edge != null; edge = prev[edge.from])
			bottleNeck = min(bottleNeck, edge.remainingCapacity());

		// Retrace augmented path and update flow values.
		for(Edge edge = prev[bigT]; edge != null; edge = prev[edge.from])
			edge.augment(bottleNeck);

		// Return bottleneck flow
		return bottleNeck;
	}


	/* Example */
	//"C:\\Users\\JERLocal\\eclipse-workspace\\testing\\src\\testing\\alex.txt"
	public static void main(String[] args) throws IOException, InterruptedException {
		long startTimeBFS = 0;//System.nanoTime();
		long endTimeBFS = 0;//System.nanoTime();
		long startTimeDFS = 0;//System.nanoTime();
		long endTimeDFS = 0;//System.nanoTime();
		testSmallFlowGraph();

	}

	// Testing graph from:
	// http://crypto.cs.mcgill.ca/~crepeau/COMP251/KeyNoteSlides/07demo-maxflowCS-C.pdf
	private static void testSmallFlowGraph() {

		//List<Edge>[] originalGraph = ;
		int n = 9;
		int s = 0;
		int t = 8;

		EdmondsKarpAdjacencyList solver;
		solver = new EdmondsKarpAdjacencyList(n, s, t);


		solver.addEdge(0, 1, 1);
		solver.addEdge(0, 2, 1);
		solver.addEdge(0, 3, 1);

		solver.addEdge(1, 4, 1);
		solver.addEdge(2, 4, 1);
		solver.addEdge(3, 4, 1);
		
		solver.addEdge(4, 5, 1);
		solver.addEdge(4, 6, 1);
		solver.addEdge(4, 7, 1);
		
		solver.addEdge(5, 8, 1);
		solver.addEdge(6, 8, 1);
		solver.addEdge(7, 8, 1);
		
		/*
		int n = 5;
		int s = 1;
		int t = 5;

		EdmondsKarpAdjacencyList solver;
		solver = new EdmondsKarpAdjacencyList(n, s, t);

		// Source edges
		solver.addEdge(1, 2, 1);
		solver.addEdge(1, 3, 1);
		solver.addEdge(1, 4, 1);
		
		solver.deepGraphBackup(solver.bigGraph);*/

		System.out.println("Max flow is: " + solver.getMaxFlow(true));
		System.out.println("Nodes Left of cut:");
		boolean[] x = solver.getMinCut(true); 
		for(int i = 0; i < x.length; i++){
			if(x[i])
				System.out.println(i);
		}

		
//		System.out.println("");
//		List<Edge>[] graph = solver.getGraph(false);
//		for(int i = 0; i < graph.length; i++){
//			for(Edge e : graph[i]){
//				if(e.isMinCut() && !e.isResidual())
//					System.out.println(e.toString(s, t));
//			}
//		}
		
		System.out.println("Cut Nodes in original graph:");
		for(int node : originalCutNodes) {
			System.out.println(node);			
		}
		System.out.println("Cut Edges in big graph:");
		for(Edge e : cutEdges) {
			System.out.println(e.toString(s, t));			
		}
	}

	/*
	 * Build graph from file, first line is num of vertices, source vertex num, target vertex num
	 * all other lines are edges: from vertex num, to vertex num, capacity of the edge
	 */
	private static EdmondsKarpAdjacencyList buildGraphfromFile(String filePath) throws IOException {
		EdmondsKarpAdjacencyList solver;
		LineNumberReader lnr = new LineNumberReader(new FileReader(filePath));
		lnr.setLineNumber(1);
		StreamTokenizer stok = new StreamTokenizer(lnr);
		stok.parseNumbers();
		stok.eolIsSignificant(true);
		stok.nextToken();
		int vertices = 0, source = 0, target = 0;
		//parse first line
		while (stok.ttype != StreamTokenizer.TT_EOL) {
			if (stok.ttype == StreamTokenizer.TT_NUMBER)
				vertices = (int)stok.nval;
			stok.nextToken();
			if (stok.ttype == StreamTokenizer.TT_NUMBER)
				source = (int)stok.nval;
			stok.nextToken();
			if (stok.ttype == StreamTokenizer.TT_NUMBER)
				target = (int)stok.nval;
			stok.nextToken();
		}
		//System.out.println("new graph, edges: " + vertices + " source: " + source + " target: " + target);
		solver = new EdmondsKarpAdjacencyList(vertices, source, target);
		stok.nextToken();
		int from = 0, to = 0, capacity = 0;
		//parse all other lines, each line is edge
		while (stok.ttype != StreamTokenizer.TT_EOF) {
			int lineNum = lnr.getLineNumber();
			double sum = 0;
			while (stok.ttype != StreamTokenizer.TT_EOL) {
				if (stok.ttype == StreamTokenizer.TT_NUMBER)
					from = (int)stok.nval;
				stok.nextToken();
				if (stok.ttype == StreamTokenizer.TT_NUMBER)
					to = (int)stok.nval;
				stok.nextToken();
				if (stok.ttype == StreamTokenizer.TT_NUMBER)
					capacity = (int)stok.nval;
				stok.nextToken();
			}
			//System.out.println("edge in line: " + lineNum + " from: " + from + " to: " + to + " capacity: " + capacity);
			if(from == source || to == target)
				capacity = Integer.MAX_VALUE / 2;
			solver.addEdge(from, to, capacity);
			stok.nextToken();
		}
		lnr.close();
		return solver;

	}

	/*
	 * Create random graph of given size and represent it in a file
	 */
	private static Boolean fillFile(String filePath, int vertices) throws IOException, InterruptedException  {
		Path fileToDeletePath = Paths.get(filePath);
		Files.deleteIfExists(fileToDeletePath);
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(filePath, true));
		}	
		catch( java.io.FileNotFoundException ex) {
			System.out.println("Access denied will retry in 3 sec");
			Thread.sleep(3000);
			writer = new BufferedWriter(new FileWriter(filePath, true));
		}


		String line = "";
		line = line.concat("" + vertices + " "); //how many vertices
		line = line.concat("" + 0 + " "); //source
		line = line.concat("" + (vertices - 1) + "\n"); //target
		writer.write(line);


		Random rand = new Random();
		List < Integer > usedVertices = new ArrayList <  > ();
		int numOfEdges = 0;
		int toVertex = 0;
		int rightMostVertex = 0;
		int leftMostVertex = 0;
		int capacity = 0;
		int possibleEdges = 20;//vertices / ; //how many edges max go out from one vertex 
		for(int fromVertex = 0; fromVertex < vertices - 1; fromVertex++) { //go over all edges except target
			//decide how many edges will go from this vertex
			numOfEdges = possibleEdges;//rand.nextInt(possibleEdges) + 1; //no more than quarter of max possible edges and at least 1 edge 
			usedVertices.add(fromVertex); //add current fromVertex to list of used vertices to avoid self pointing
			//limits for edges going from this vertex, go more right than left?
			rightMostVertex = Math.min((fromVertex + (int)(possibleEdges)), (vertices - 1));
			leftMostVertex = Math.max(fromVertex - (int)(possibleEdges), 0);
			//add to file all edges going from current fromVertex
			for(int edge = 0; edge < numOfEdges; edge++) { 
				do {//choose the destination of edge, can be any vertex between leftmost and rightmost
					toVertex = rand.nextInt(rightMostVertex - leftMostVertex + 1) + leftMostVertex; 
				}while(usedVertices.contains(toVertex));
				usedVertices.add(toVertex); //add the toVertex to the list of used vertices, to avoid multiple edges between same nodes
				line = "";
				line = line.concat("" + fromVertex + " ");
				line = line.concat("" + toVertex + " ");
				capacity =  1;
				line = line.concat("" + capacity + "\n"); //capacity is always 1?
				writer.write(line);

				//System.out.println(line);			

			}
			//last quarter of vertex always connect to target
			//			if((fromVertex >= vertices - possibleEdges) && !usedVertices.contains(vertices - 1)) {
			//				System.out.println("last quarter");
			//				line = "";
			//				line = line.concat("" + fromVertex + " ");
			//				line = line.concat("" + (vertices - 1) + " ");
			//				capacity = rand.nextInt(10) + 1;
			//				line = line.concat("" + capacity + "\n"); //capacity is always 1?
			//				System.out.println("last quarter line " + line);
			//				writer.write(line);
			//			}
			usedVertices.clear();//clean the list of used vertices before moving to next fromVertex
		}
		line = "\n";
		writer.write(line);
		writer.close();
		return true;
		//System.out.println("done filling file");
	}
}