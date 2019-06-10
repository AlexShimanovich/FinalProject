// Kosaraju's algorithm for finding  SCC in a graph
import java.io.*; 
import java.util.*; 
import java.util.LinkedList; 


class GreedGraph 
{ 
	private int V; // No. of vertices in graph 
	private LinkedList<Integer> adj[]; //Adjacency List Graph representation
	int rootArray[]; //track the origin of each node (same origin == same SCC)
	int sizeSccArray[]; //measure size of SCC for each origin vertex, left zero if this vertex is not origin

	//Constructor 
	GreedGraph(int v){ 
		V = v; 
		adj = new LinkedList[v]; 		
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
			if(!visited[n]) 
				fillDfsOrder(n, visited, stack); 
		} 

		// All vertices reachable from v are processed by now, 
		// push v to Stack 
		stack.push(new Integer(v)); 
	} 
	
	// A recursive function to print DFS starting from v 
	private void DFSUtil(int v, boolean visited[], int origin){ 
		// Mark the current node as visited and print it 
		visited[v] = true;
		rootArray[v] = origin;
		sizeSccArray[origin]++;
		System.out.print(v + " "); 

		int n; 

		// Recursion for all the vertices adjacent to this vertex 
		Iterator<Integer> i = adj[v].iterator(); 
		while (i.hasNext()) 
		{ 
			n = i.next(); 
			if (!visited[n]) 
				DFSUtil(n, visited, origin); 
		} 
	} 

	// The main function that finds and prints all strongly 
	// connected components 
	private void printSCCs() 
	{ 
		Stack stack = new Stack(); 

		// Mark all the vertices as not visited (For first DFS) 
		boolean visited[] = new boolean[V]; 
		for(int i = 0; i < V; i++) 
			visited[i] = false; 

		// Fill vertices in stack according to their finishing times 
		for (int i = 0; i < V; i++) 
			if (!visited[i]) 
				fillDfsOrder(i, visited, stack); 

		// Create a transposed graph 
		GreedGraph gr = getTranspose(); 

		// Mark all the vertices as not visited (For second DFS) 
		for (int i = 0; i < V; i++) 
			visited[i] = false; 

		
		// Now process all vertices in order defined by Stack 
		while (!stack.empty()){ 
			// Pop a vertex from stack 
			int v = (int)stack.pop(); 
			//root[v]++;
			// Print Strongly connected component of the popped vertex 
			if (!visited[v]) 
			{ 
				gr.DFSUtil(v, visited, v); //running on transposed graph
				System.out.println(); 
			} 
		}
		System.out.println("Node Origins array: " + Arrays.toString(gr.rootArray));
		System.out.println("SCC Size array: " + Arrays.toString(gr.sizeSccArray));
		Arrays.sort(gr.sizeSccArray);
		System.out.println("Biggest SCC size is " + gr.sizeSccArray[gr.sizeSccArray.length - 1]);
	} 

	// Driver method 
	public static void main(String args[]) 
	{ 
		// Create a graph given in the above diagram 
		GreedGraph g = new GreedGraph(5); 
		g.addEdge(1, 0); 
		g.addEdge(0, 2); 
		g.addEdge(2, 1); 
		g.addEdge(0, 3); 
		g.addEdge(3, 4); 

		System.out.println("Those are the SCC of our graph:"); 
		g.printSCCs(); 
	} 
} 
