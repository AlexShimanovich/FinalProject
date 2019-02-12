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

//////////////////////////////////Original BFS///////////////////////////////////////	
	
	
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
	
//////////////////////////////////Big BFS///////////////////////////////////////	
	
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
				//convert from bigLeftOfCut to leftOfCut
				if(i % 2 == 0)
					leftOfCut[i / 2] = true;
				
				
				//update cut edges, if edge goes to visited node its not in cut
				ArrayList<Edge> edgesToRemove =  new ArrayList<Edge>();
				for(Edge edge : cutEdges) {
					if(edge.to == i)
						edgesToRemove.add(edge);
				}
				cutEdges.removeAll(edgesToRemove);
			}
		}
		
		//were done, convert from big edges to original graph nodes
		for(Edge edge : cutEdges) {
			originalCutNodes.add(edge.from / 2);
		}
				
	}
	
	private long bigBfs(long CurrentFlow) {
		Edge[] prev = new Edge[bigN];
		cutEdges.clear();
		// The queue can be optimized to use a faster queue
		Queue<Integer> q = new ArrayDeque<>(bigN);
		bigVisit(bigS);
		q.offer(bigS);
		
		//System.out.println("START bfsssssssssss, current flow: " + CurrentFlow );
		// Perform BFS from source to sink
		while(!q.isEmpty()) {
			int node = q.poll();
			if (node == bigT) break;
			//System.out.println("node " + node );
			for (Edge edge : bigGraph[node]) {
				long cap = edge.remainingCapacity();
				if (cap > 0 && !bigVisited(edge.to)) {
					bigVisit(edge.to);
					prev[edge.to] = edge;
					q.offer(edge.to);
				}
				//System.out.println("edge " + edge.toString(edge.from, edge.to) );
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
	
	//calculate current flow through the cut
	private long floww() {
		long flow = 0;
		for (Edge e : cutEdges)
			flow += e.capacity;
		return flow;	
	}


//////////////////////////////////MAIN///////////////////////////////////////	
	
	/* Example */
	//"C:\\Users\\JERLocal\\eclipse-workspace\\testing\\src\\testing\\alex.txt"
	public static void main(String[] args) throws IOException, InterruptedException {
		long startTimeBFS = 0;//System.nanoTime();
		long endTimeBFS = 0;//System.nanoTime();
		long startTimeDFS = 0;//System.nanoTime();
		long endTimeDFS = 0;//System.nanoTime();
		testSmallFlowGraph();
	}


	private static void testSmallFlowGraph() throws IOException {
		
		List<Edge>[] originalGraph = graphfromFile("C:\\Users\\JERLocal\\eclipse-workspace\\testing\\src\\testing\\choke.txt");
		
		//print graph
//		for (int i = 0; i < originalGraph.length; i++) {
//			System.out.println("graph from text node i = " + i);
//			for (Edge e : originalGraph[i]) {
//				System.out.println(e.toString(e.from, e.to));				
//			}					
//		}
		List<Integer> W = new ArrayList<Integer>();
		List<Integer> A = new ArrayList<Integer>();
		List<Integer> B = new ArrayList<Integer>();
		
		for (int wSize = 4; wSize < 16; wSize *= 2) {
			W = getRandomWFromGraph(originalGraph.length, wSize);
			A = new ArrayList<Integer>();
			B = new ArrayList<Integer>();
			System.out.println("//////////////////////////////////// W size is: " + wSize + "////////////////////////////////");
			System.out.println("W is:");
			for(int x : W) {
				System.out.print(x + ",");
			}
			System.out.println("");
			int t = (int) Math.min(Math.pow(2, wSize), 256);
			
			for(int repeatSameW = 0; repeatSameW < t; repeatSameW ++){	
				System.out.println("////////////////// W size is: " + wSize + ", t is: " + repeatSameW + " out of " + (t-1) );
				fillAandB(W, A, B);
				System.out.println("A is:");
				for(int x : A) {
					System.out.print(x + ",");
				}
				System.out.println("");
				System.out.println("B is:");
				for(int x : B) {
					System.out.print(x + ",");
				}
				System.out.println("");			
				EdmondsKarpAdjacencyList solver = solverFromGraph(originalGraph, W, A, B);
				long maxFlow = solver.getMaxFlow(true);
				System.out.println("Max flow is: " + maxFlow);
				if(maxFlow == 0)
					continue; //no flow, continue to another A and B

				System.out.println("Nodes Left of cut in BigGraph:");
				boolean[] x = solver.getMinCut(true); 
				for(int i = 0; i < x.length - 2; i++){
					if(x[i])
						System.out.print(i + ",");
				}
				System.out.println("");
				
				System.out.println("Nodes LEFT of cut in W:");
				for(int i = 0; i < solver.leftOfCut.length - 2; i++){  // -2 because we dont want s,t
					if(W.contains(i) && solver.leftOfCut[i] && !solver.originalCutNodes.contains(i))
						System.out.print(i + ",");
				}
				System.out.println("");
				
				System.out.println("Nodes RIGHT of cut in W:");
				for(int i = 0; i < solver.leftOfCut.length - 2; i++){ // -2 because we dont want s,t
					if(W.contains(i) && !solver.leftOfCut[i] && !solver.originalCutNodes.contains(i) )
						System.out.print(i + ",");
				}		
				System.out.println("");
				
				System.out.println("Cut Nodes in W:");
				for(int node : solver.originalCutNodes) {
					System.out.println(node);			
				}
				System.out.println("Cut Edges in big graph:");
				for(Edge e : solver.cutEdges) {
					System.out.println(e.toString(e.from, e.to));			
				}	
			}


			
		}

		
		//random
//		List<Integer> W = new ArrayList<Integer>(Arrays.asList(0,1,2,3,4,5,6));
//		List<Integer> A = new ArrayList<Integer>(Arrays.asList(0,1,2));
//		List<Integer> B = new ArrayList<Integer>(Arrays.asList(3,4,5,6));
//		
//		EdmondsKarpAdjacencyList solver = solverFromGraph(originalGraph, W, A, B);
//		
//		System.out.println("Max flow is: " + solver.getMaxFlow(true));
//		System.out.println("Nodes Left of cut in BigGraph:");
//		boolean[] x = solver.getMinCut(true); 
//		for(int i = 0; i < x.length - 2; i++){
//			if(x[i])
//				System.out.print(i + ",");
//		}
//		System.out.println("");
//		
//		System.out.println("Nodes LEFT of cut in W:");
//		for(int i = 0; i < solver.leftOfCut.length - 2; i++){  // -2 because we dont want s,t
//			if(W.contains(i) && solver.leftOfCut[i] && !originalCutNodes.contains(i))
//				System.out.print(i + ",");
//		}
//		System.out.println("");
//		
//		System.out.println("Nodes RIGHT of cut in W:");
//		for(int i = 0; i < solver.leftOfCut.length - 2; i++){ // -2 because we dont want s,t
//			if(W.contains(i) && !solver.leftOfCut[i] && !originalCutNodes.contains(i) )
//				System.out.print(i + ",");
//		}		
//		System.out.println("");
//		
//		System.out.println("Cut Nodes in W:");
//		for(int node : originalCutNodes) {
//			System.out.println(node);			
//		}
//		System.out.println("Cut Edges in big graph:");
//		for(Edge e : cutEdges) {
//			System.out.println(e.toString(e.from, e.to));			
//		}

	}

	private static List<Integer> getRandomWFromGraph(int nodesRange, int wSize){
		List<Integer> W = new ArrayList<Integer>();
		Random r = new Random();
		int randomVertex;
		while(W.size() < wSize) {
			randomVertex = r.nextInt((nodesRange)); //ints from [0,nodesRange)
			if( !W.contains(randomVertex))
				W.add(randomVertex);			
		}
		return W;		
	}
	
	//fill A and B randomly with nodes from W
	//if A or B reach 0.75 of W size , put rest of nodes in B or A accordingly
	private static void fillAandB(List<Integer> W, List<Integer> A, List<Integer> B) {
		A.clear();
		B.clear();
		for(int node : W) {
			if(A.size() < (0.75 * W.size()) && B.size() < (0.75 * W.size())) {
				double x = Math.random();

				if(x > 0.5)
					A.add(node);
				else
					B.add(node);
				continue;
			}

			//A or B reached 0.75 * W			
			if(B.size() >= (0.75 * W.size())) {
				A.add(node);
				continue;
			}
			else {
				B.add(node);	
			}

		}
				
		
	}
	
	/*
	 * create solver from given graph containing only nodes from A connected to source and nodes From b connected to target
	 */
	private static EdmondsKarpAdjacencyList solverFromGraph(List<Edge>[] originalGraph, List<Integer> W, List<Integer> A, List<Integer> B) {
		int maxNodeInW = W.stream()
			      .mapToInt(v -> v)
			      .max().orElseThrow(NoSuchElementException::new);
		//System.out.println("max node innnnnnn" + maxNodeInW);
		//solver from originalGraph + random, solver size w+2, source at location w, target at location w+1
		EdmondsKarpAdjacencyList solver;
		int verticesSolver = maxNodeInW + 3;
		int s = maxNodeInW + 1;
		int t = maxNodeInW + 2;
		solver = new EdmondsKarpAdjacencyList(verticesSolver,  s, t);
		//copy existing edges but only those  that both vertices are in w
		for (int i : W) {
			//System.out.println("building solver i = " + i);
			for (Edge e : originalGraph[i]) {
				if(W.contains(e.to)) {
					solver.addEdge(e.from, e.to, e.capacity);
					//System.out.println(e.toString(e.from, e.to));
				}
			}					
		}
		//add edges from s to A with big capacity
		for(int a : A) {
			solver.addEdge(s, a, verticesSolver);
		}
		//add edges from B to t with big capacity
		for(int b : B) {
			solver.addEdge(b, t, verticesSolver);
		}
		
		return solver;
		
	}

//////////////////////////////////File methods///////////////////////////////////////
	
	
	/*
	 * Build solver from file, first line is num of vertices, source vertex num, target vertex num
	 * all other lines are edges: from vertex num, to vertex num, capacity of the edge
	 */
	private static EdmondsKarpAdjacencyList buildSolverfromFile(String filePath) throws IOException {
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
		//System.out.println("new graph, vertices: " + vertices + " source: " + source + " target: " + target);
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
	 * Build graph from file, first line is num of vertices
	 * all other lines are edges: from vertex num, to vertex num, capacity of the edge
	 * last lien empty
	 */
	private static List<Edge>[] graphfromFile(String filePath) throws IOException {
		LineNumberReader lnr = new LineNumberReader(new FileReader(filePath));
		lnr.setLineNumber(1);
		StreamTokenizer stok = new StreamTokenizer(lnr);
		stok.parseNumbers();
		stok.eolIsSignificant(true);
		stok.nextToken();
		int vertices = 0;
		//parse first line
		while (stok.ttype != StreamTokenizer.TT_EOL) {
			if (stok.ttype == StreamTokenizer.TT_NUMBER) {
				vertices = (int)stok.nval;
				stok.nextToken();
				break;
			}
			//file ERROR

		}
//		System.out.println(vertices);
		List<Edge>[] originalGraph = new List[vertices];
		for (int i = 0; i < originalGraph.length; i++)
			originalGraph[i] = new ArrayList<Edge>();
		stok.nextToken();
		int from = 0, to = 0, capacity = 0;
		//parse all other lines, each line is edge
		while (stok.ttype != StreamTokenizer.TT_EOF) {
			int lineNum = lnr.getLineNumber();
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
//			System.out.println("edge in line: " + lineNum + " from: " + from + " to: " + to + " capacity: " + capacity);
			originalGraph[from].add(new Edge(from, to, capacity));
			stok.nextToken();
		}
		lnr.close();
		return originalGraph;

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