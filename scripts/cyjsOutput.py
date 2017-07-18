from collections import Counter
from collections import defaultdict

from porestat.utils.DataFrame import DataFrame, ExportTYPE

interactions = DataFrame.parseFromFile("/home/users/joppich/ownCloud/data/chemokines_sfb/chemokine_interactions.tsv" )

#interactions.export("/home/users/joppich/ownCloud/data/chemokines_sfb/chemokine_interactions.html", ExportTYPE.HTML)

uniqueEdges = set()

chemList = [
    'CXCR2',
'CCL9',
'CXCL5',
'CXCL1',
'CXCL13',
'CXCL7',
'CCL2',
'CXCL9',
'CCL3',
'CXCL10',
'CCL22',
'CCR5',
'CCR7',
'CCL7',
'CCL4',
'CXCR4',
'CX3CL1',
'CXCL12',
'CXCL14',

]

cnt = 0
usedSources = set()

for row in interactions:

    if row.data[0] == 'tarbase_v7':
        continue

    if row.data[4] > 0.95 or row.data[0] == 'coocurrence':

        if row.data[0] == 'targetscan_v7' or row.data[0] == 'tarbase_v7':
            if row.data[4] < 0.99:
                continue

        mirname = "-".join(row.data[6].upper().split("-")[1:])
        genename = row.data[5].upper()

        if not genename in chemList:
            continue

        usedSources.add(row.data[0])

        edge = (genename, mirname)

        uniqueEdges.add( edge )


mirsByGene = defaultdict(set)

for edge in uniqueEdges:

    print( edge )
    mirsByGene[edge[0]].add(edge[1])

print(len(uniqueEdges))

nodeID = 0
allNodes = {}

for gene in mirsByGene:

    if not gene in allNodes:
        allNodes[gene] = 'n'+str(nodeID)
        nodeID += 1

    for mirna in mirsByGene [gene]:
        if mirna in allNodes:
            continue

        allNodes[mirna] = 'n'+str(nodeID)
        nodeID += 1


elementsStr = ""

elementsStr += "elements: {\n nodes: ["

for node in allNodes:

    color = '#11479e'

    if node.startswith("MIR") or node.startswith("LET"):
        color='#bf3d3d'

    elementsStr += "{ data: { id: '" + allNodes[node] + "', name: '"+node+"', color: '"+color+"'}},"


elementsStr += "], \n edges: ["

for edge in uniqueEdges:
    sourceid = allNodes[edge[0]]
    targetid = allNodes[edge[1]]

    elementsStr += "{ data: { source: '"+sourceid+"', target: '"+targetid+"'}},"

elementsStr += "]\n}"


HTMLtop = """
<!DOCTYPE>
<!-- This code is for demonstration purposes only.  You should not hotlink to Github, Rawgit, or files from the Cytoscape.js documentation in your production apps. -->
<html>
  <head>
    <title>cose demo</title>
    <meta name="viewport" content="width=device-width, user-scalable=no, initial-scale=1, maximum-scale=1">
    <!-- For loading external data files -->
    <script src="http://www.trauen-sich.net/js/cytoscape.js"></script>

	<style>

body {
  font-family: helvetica;
  font-size: 14px;
}

#cy {
  width: 100%;
  height: 100%;
  position: absolute;
  left: 0;
  top: 0;
  z-index: 999;
}

h1 {
  opacity: 0.5;
  font-size: 1em;
}

</style>

  </head>
  <body>
    <h1>Chemokines</h1>
    <div id="cy"></div>
    <!-- Load appplication code at the end to ensure DOM is loaded -->
    <script>
console.log("hello");

    var cy = window.cy = cytoscape({
      container: document.getElementById('cy'),

boxSelectionEnabled: false,
  autounselectify: true,

  layout: {
  name: 'cose',
  // Called on `layoutready`
  ready: function(){},
  // Called on `layoutstop`
  stop: function(){},
  // Whether to animate while running the layout
  animate: true,
  // The layout animates only after this many milliseconds
  // (prevents flashing on fast runs)
  animationThreshold: 250,
  // Number of iterations between consecutive screen positions update
  // (0 -> only updated on the end)
  refresh: 20,
  // Whether to fit the network view after when done
  fit: true,
  // Padding on fit
  padding: 30,
  // Constrain layout bounds; { x1, y1, x2, y2 } or { x1, y1, w, h }
  boundingBox: undefined,
  // Excludes the label when calculating node bounding boxes for the layout algorithm
  nodeDimensionsIncludeLabels: false,
  // Randomize the initial positions of the nodes (true) or use existing positions (false)
  randomize: false,
  // Extra spacing between components in non-compound graphs
  componentSpacing: 100,
  // Node repulsion (non overlapping) multiplier
  nodeRepulsion: function( node ){ return 400000; },
  // Node repulsion (overlapping) multiplier
  nodeOverlap: 10,
  // Ideal edge (non nested) length
  idealEdgeLength: function( edge ){ return 40; },
  // Divisor to compute edge forces
  edgeElasticity: function( edge ){ return 100; },
  // Nesting factor (multiplier) to compute ideal edge length for nested edges
  nestingFactor: 5,
  // Gravity force (constant)
  gravity: 80,
  // Maximum number of iterations to perform
  numIter: 1000,
  // Initial temperature (maximum node displacement)
  initialTemp: 200,
  // Cooling factor (how the temperature is reduced between consecutive iterations
  coolingFactor: 0.95,
  // Lower temperature threshold (below this point the layout will end)
  minTemp: 1.0,
  // Pass a reference to weaver to use threads for calculations
  weaver: false
  },

  style: [
    {
      selector: 'node',
      style: {
        'content': 'data(name)',
        'text-opacity': 0.5,
        'text-valign': 'center',
        'text-halign': 'right',
        'background-color': 'data(color)'
      }
    },

    {
      selector: 'edge',
      style: {
        'width': 4,
        'target-arrow-shape': 'triangle',
        'line-color': '#9dbaea',
        'target-arrow-color': '#9dbaea'
      }
    }
  ],
"""

HTMLbottom = """


    });

</script>
  </body>
</html>
"""

with open("cyjs_plot.html", 'w') as outHTML:

    outHTML.write( HTMLtop + elementsStr + HTMLbottom );