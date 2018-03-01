var TreeModel = require('tree-model');


var tree = new TreeModel();


function Node(value) {

    this.value = value;
    this.children = [];
    this.parent = null;

    this.setParentNode = function(node) {
        this.parent = node;
    }

    this.getParentNode = function() {
        return this.parent;
    }

    this.addChild = function(node) {
        node.setParentNode(this);
        this.children[this.children.length] = node;
    }

    this.getChildren = function() {
        return this.children;
    }

    this.removeChildren = function() {
        this.children = [];
    }
}

var root = new Node('root');
root.addChild(new Node('child 0'));
root.addChild(new Node('child 1'));
var children = root.getChildren();
for(var i = 0; i < children.length; i++) {
    for(var j = 0; j < 5; j++) {
        children[i].addChild(new Node('second level child ' + j));
    }
}
console.log(root);
children[0].removeChildren();
console.log(root);
console.log(root.getParentNode());
console.log(children[1].getParentNode());

var spheres = {
	s1:[5],
	s2:[4,6,7],
	s3:[[3],[1],[8,9]],
	s4:[[0,2],[2],[],[10,15]],
	s5:[[],[1,16],[11,12],[14]]
};


/*
root = tree.parse({
    id: 5,
    children: [
        {
            id: 4,
            children: [{id: 3}]
        },
        {
            id: 6,
            children: [{id: 1}]
        },
        {
            id: 7,            
            children: [{id: 8},{id: 9}]

        }
    ]
});
*/

/*
var t = 0;
for (obj in spheres){
	  t=t+1;

	  if (t==1) root = tree.parse({id: t-1, name: obj, bondtype: 1.5, symbol: 'C', rank: 1, score: 0});
	  for (j=0;j<obj.length;j++) {
		  var subnode = tree.parse({id: t-1, name: obj, bondtype: 1.5, symbol: 'C', rank: 1, score: 0}); // Parse the new node
		  console.log(subnode);		  
		  root.addChild(subnode); 		// Add it to the parent
		}
}
*/
//console.log(root.getPath());