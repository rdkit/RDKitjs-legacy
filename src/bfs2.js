/* A Queue object for queue-like functionality over JavaScript arrays. */
var Queue = function() {
    this.items = [];
};
Queue.prototype.enqueue = function(obj) {
    this.items.push(obj);
};
Queue.prototype.dequeue = function() {
    return this.items.shift();
};
Queue.prototype.isEmpty = function() {
    return this.items.length === 0;
};

/*
 * Performs a breadth-first search on a graph
 * @param {array} graph - Graph, represented as adjacency lists.
 * @param {number} source - The index of the source vertex.
 * @returns {array} Array of objects describing each vertex, like
 *     [{distance: _, predecessor: _ }]
 */
var doBFS = function(graph, source) {
    var bfsInfo = [];

    for (var i = 0; i < graph.length; i++) {
	    bfsInfo[i] = {
	        distance: null,
	        predecessor: null };
    }

    bfsInfo[source].distance = 0;

    var queue = new Queue();
    queue.enqueue(source);

    // Traverse the graph
    
    // As long as the queue is not empty:
    //  Repeatedly dequeue a vertex u from the queue.
    //  
    //  For each neighbor v of u that has not been visited:
    //     Set distance to 1 greater than u's distance
    //     Set predecessor to u
    //     Enqueue v
    //
    //  Hint:
    //  use graph to get the neighbors,
    //  use bfsInfo for distances and predecessors 
    
    return bfsInfo;
};


var adjList = [
    [1],
    [0, 4, 5],
    [3, 4, 5],
    [2, 6],
    [1, 2],
    [1, 2, 6],
    [3, 5],
    []
    ];
var bfsInfo = doBFS(adjList, 3);
for (var i = 0; i < adjList.length; i++) {
    println("vertex " + i + ": distance = " + bfsInfo[i].distance + ", predecessor = " + bfsInfo[i].predecessor);
}

/*
Program.assertEqual(bfsInfo[0], {distance: 4, predecessor: 1});
Program.assertEqual(bfsInfo[1], {distance: 3, predecessor: 4});
Program.assertEqual(bfsInfo[2], {distance: 1, predecessor: 3});
Program.assertEqual(bfsInfo[3], {distance: 0, predecessor: null});
Program.assertEqual(bfsInfo[4], {distance: 2, predecessor: 2});
Program.assertEqual(bfsInfo[5], {distance: 2, predecessor: 2});
Program.assertEqual(bfsInfo[6], {distance: 1, predecessor: 3});
Program.assertEqual(bfsInfo[7], {distance: null, predecessor: null});
*/