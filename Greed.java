// Kosaraju's (DFS based) algorithm for finding  SCC in a graph
//Alex Shimanovich

import java.io.*; 
import java.util.*; 
import java.util.LinkedList; 


class GreedGraph 
{ 

	private static class Result{
		private int SCC_Size;
		private LinkedList<Integer> indexesOfSeparator;		

		public Result(int SCC_Size, LinkedList<Integer> indexesOfSeparator) {
			this.SCC_Size = SCC_Size;
			this.indexesOfSeparator = new LinkedList<Integer>(indexesOfSeparator); 
		}

		public int getSccSize() {
			return SCC_Size;
		}

		public LinkedList<Integer> getIndexesOfSeparator() {
			return indexesOfSeparator;
		}

		public int getSeparatorSize() {
			return indexesOfSeparator.size();	
		}

		public String toString() {
			String res = "SCC size: " + SCC_Size + "\n Separator size: " + indexesOfSeparator.size();  
			return res;
		}

	}

	private int V; // No. of vertices in graph 
	private LinkedList<Integer> adj[]; //Adjacency List Graph representation
	private int rootArray[]; //track the origin of each node (same origin == same SCC)
	private int sizeSccArray[]; //measure size of SCC for each origin vertex, left zero if this vertex is not origin
	private int separatorArray[]; //nodes from graph that are part of separator
	private LinkedList<Integer> indexesOfSeparator; //list of indexes (nodes) that are in separator
	//token to be used to know if node is part of separator, 
	//instead of changing 0 to 1 and vice-versa in all the separatorArray[] I just do separatorToken++
	private int separatorToken = 0; 

	//Constructor 
	GreedGraph(int v){ 
		V = v; 
		adj = new LinkedList[V]; 
		separatorArray = new int[V];
		indexesOfSeparator = new LinkedList();
		for (int i = 0; i < v; ++i) 
			adj[i] = new LinkedList(); 
	} 

	/**
	 * Function to add an edge into the graph 
	 * @param v source vertex
	 * @param w destination vertex
	 */
	private void addEdge(int v, int w){ 
		adj[v].add(w); 
	} 

	// Function that returns transpose of this graph 
	private GreedGraph getTranspose() 
	{ 
		GreedGraph transposeGraph = new GreedGraph(V); 
		for (int v = 0; v < V; v++) 
		{ 
			// Recursion for all the vertices adjacent to this vertex 
			Iterator<Integer> i = adj[v].listIterator(); 
			while(i.hasNext()) 
				transposeGraph.adj[i.next()].add(v); 
		} 
		transposeGraph.rootArray = new int[V];
		transposeGraph.sizeSccArray = new int[V];
		transposeGraph.separatorArray = separatorArray;
		transposeGraph.indexesOfSeparator = indexesOfSeparator;
		transposeGraph.separatorToken = separatorToken;
		return transposeGraph; 
	} 

	// Fill vertices in stack according to their DFS finishing times 
	private void fillDfsOrder(int v, boolean visited[], Stack stack) 
	{ 
		// Mark the current node as visited
		visited[v] = true; 

		// Recursion for all the vertices adjacent to this vertex 
		Iterator<Integer> i = adj[v].iterator(); 
		while (i.hasNext()) 
		{ 
			int n = i.next(); 
			if(!visited[n] /*&& separatorArray[n] != separatorToken*/) 
				fillDfsOrder(n, visited, stack); 
		} 

		// All vertices reachable from vertex v are processed by now, 
		// push v to Stack 
		stack.push(new Integer(v)); 
	} 

	// A recursive function to print DFS starting from v 
	private void DFSUtil(int v, boolean visited[], int origin, boolean print){ 
		// Mark the current node as visited and print it 
		visited[v] = true;
		rootArray[v] = origin;
		sizeSccArray[origin]++;
		if(print)
			System.out.print(v + " "); 

		// Recursion for all the vertices adjacent to this vertex 
        for(Integer node: adj[v]){ 
        	//System.out.print("looking at node " + node);
			if (!visited[node] && separatorArray[node] != separatorToken) 
				DFSUtil(node, visited, origin, print); 
        } 

	} 

	// The main function that finds and prints all strongly 
	// connected components in Graph
	// return the size of biggest SCC in the graph
	private int findSCCs(boolean print) 
	{ 
		Stack stack = new Stack(); 

		// Mark all the vertices as not visited (For first DFS) 
		boolean visited[] = new boolean[V]; 
		for(int i = 0; i < V; i++) 
			visited[i] = false; 

		// Fill vertices in stack according to their finishing times 
		for (int i = 0; i < V; i++) 
			if (!visited[i] /*&& separatorArray[i] != separatorToken*/) 
				fillDfsOrder(i, visited, stack); 

		// Create a transposed graph 
		GreedGraph gr = getTranspose(); 

		// Mark all the vertices as not visited (For second DFS) 
		for (int i = 0; i < V; i++) 
			visited[i] = false; 
/*
        for(int v = 0; v < gr.V; v++) 
        { 
            System.out.println("Adjacency list of vertex "+ v); 
            System.out.print("head"); 
            for(Integer pCrawl: gr.adj[v]){ 
                System.out.print(" -> "+pCrawl); 
            } 
            System.out.println("\n"); 
        } */

		// Now process all vertices in transpose graph in order defined by Stack 
		while (!stack.empty()){ 
			// Pop a vertex from stack 
			int v = (int)stack.pop(); 
			//root[v]++;
			// Print Strongly connected component of the popped vertex 
			if (!visited[v] && separatorArray[v] != separatorToken) 
			{ 
				gr.DFSUtil(v, visited, v, print); //running on transposed graph
				if(print)
					System.out.println(); //new line for new SCC
			} 
		}

		Arrays.sort(gr.sizeSccArray);
		int biggestScc = gr.sizeSccArray[gr.sizeSccArray.length - 1];
//		if(print) {
//			System.out.println("Node Origins array: " + Arrays.toString(gr.rootArray));
//			System.out.println("SCC Size array: " + Arrays.toString(gr.sizeSccArray));
//			System.out.println("Biggest SCC size is " + biggestScc);
//		}
		return biggestScc;
	} 


	//get random integer in range [min, max] (inclusive)
	public int randomIntRange(int min, int max){
		int x = (int)(Math.random() * ( (max - min) + 1) ) + min;
		return x;
	}

	//randomly select n different vertices from graph to be separator, 
	//n must be smaller than V
	//we make sure that with the random separator the biggest SCC is smaller than 2/3 of V 
	public void randomSeparator(int n){
		int node;
		separatorToken++;
		indexesOfSeparator.clear();
		Random rand = new Random(); /*
		do{
			for(int i = 0; i < n; i++) {
				do {
					node = rand.nextInt(V); // random between [0,V) - every node in graph can be part of separator	
				}while(separatorArray[node] == separatorToken);	// make sure we select different node each time	
				separatorArray[node] = separatorToken; //collect all indexes of the separator
				indexesOfSeparator.add(node); //collect all indexes of the separator
			}
		}while(findSCCs(false) > (int)(0.666 * V)); //if this separator creates too big SCC look for another separator*/
		separatorArray[3] = separatorToken; //collect all indexes of the separator
		indexesOfSeparator.add(3); //collect all indexes of the separator
		separatorArray[4] = separatorToken; //collect all indexes of the separator
		indexesOfSeparator.add(4); //collect all indexes of the separator
	}

	//move one random node from the separator back to graph
	public void moveNodeBackToGraph(){
		Random rand = new Random(); 
		int randomIndex = rand.nextInt(indexesOfSeparator.size());
		int indexToMoveBack = indexesOfSeparator.get(randomIndex); //randomly select one index	
		separatorArray[indexToMoveBack]--; //remove it from the separator (it no longer equals current separator token)
		indexesOfSeparator.remove(randomIndex); //remove it from the indexes list
		System.out.println("moving back node " + indexToMoveBack);
	}

	// Test engine
	public static void main(String args[]) 
	{ 
		//test graph
		GreedGraph g = new GreedGraph(7); 
		g.addEdge(1, 0); 
		g.addEdge(0, 2); 
		g.addEdge(2, 1); 
		g.addEdge(0, 3); 
		g.addEdge(3, 4); 
		g.addEdge(4, 5);
		g.addEdge(5, 6);
		g.addEdge(6, 4);


		int iteration = 0;
		int biggestSCC = 0;
		LinkedList<Result> results = new LinkedList(); //collect the best results from each run
		Result res = null;
		while(iteration < 2) {
			int separatorSize = 3;
			g.randomSeparator(separatorSize); //create random separator in graph
			System.out.println("Start ITERATION: " + iteration); 
			System.out.println("Initial random separator is : " + Arrays.toString(g.indexesOfSeparator.toArray()));
			//each iteration runs while we can move nodes and biggest SCC is not too big 
			do{
				System.out.println("Separator now: " + Arrays.toString(g.indexesOfSeparator.toArray())); 
				g.moveNodeBackToGraph(); 
				biggestSCC = g.findSCCs(true);
				if(biggestSCC > (int)(0.666 * g.V)) {
					System.out.println("Biggest scc is too big, size: " + biggestSCC + " DEAD END");
					break;
				}
				if(biggestSCC > 1) { //collect interesting results
					res = new Result(biggestSCC, g.indexesOfSeparator);	
					results.add(res);
				}

			}while(g.indexesOfSeparator.size() > 1);

			System.out.println("End ITERATION: " + iteration);
			iteration++;
		}

		System.out.println("Lets look at some interesting results");
		if(results.size() > 0) {
			for(Result resIterator : results)
				System.out.println(resIterator);
		}


	} 
} 
