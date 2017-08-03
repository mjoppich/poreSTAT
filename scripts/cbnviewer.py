import json
from collections import Counter
from collections import defaultdict




from porestat.utils.DataFrame import DataFrame, ExportTYPE

base = "/mnt/c/ownCloud/data/bcn/"
inputFile = base + "all/CV-IPN-Endothelial cell activation1.0.sif"
inputFile = base + "manual/CV-IPN-Endothelial cell activation1.0.jgf"

allNodesLabels = {}
allEdges = []
uniqueNodes = set()

if inputFile.endswith('.sif'):
    interactions = DataFrame.parseFromFile(inputFile, ['source', 'interaction', 'target'])


    for row in interactions:

        src = row['source']
        dst = row['target']

        if src is None or dst is None:
            continue

        uniqueNodes.add(src)
        uniqueNodes.add(dst)



    for row in interactions:
        src = row['source']
        dst = row['target']
        act = row['interaction']

        if src is None or dst is None:
            continue

        allEdges.append( (src, dst, act) )

if inputFile.endswith('.jgf'):

    fileInput = open(inputFile, 'r')
    jsonGraph = json.load(fileInput)['graph']
    title = jsonGraph['label']

    jsonNodes = jsonGraph['nodes']
    jsonEdges = jsonGraph['edges']

    for nodeDetail in jsonNodes:
        nodeId = nodeDetail['id']
        nodeLabel = nodeDetail['label']

        uniqueNodes.add(nodeId)
        allNodesLabels[nodeId] = nodeLabel

    for edgeDetail in jsonEdges:
        #edgeDetail = jsonEdges[edge]

        src = edgeDetail['source']
        dst = edgeDetail['target']
        act = edgeDetail['relation']

        if src is None or dst is None:
            continue

        allEdges.append( (src, dst, act) )


allNodes = {}
nodes = list(uniqueNodes)
idx = 0
for i in range(0, len(nodes)):
    allNodes[nodes[i]] = i


elementsStr = ""
elementsStr += "elements: {\n nodes: ["

for node in allNodes:

    color = '#11479e'

    if node.startswith("MIR") or node.startswith("LET"):
        color='#bf3d3d'

    nodeLabel = node
    if node in allNodesLabels:
        nodeLabel = allNodesLabels[node]

    elementsStr += "{ data: { id: '" + str(allNodes[node]) + "', name: '"+nodeLabel+"', color: '"+color+"'}},"


elementsStr += "], \n edges: ["

for edge in allEdges:

    """
    a abundance
    g geneabundance
    m micro rna abundance
    p protein abundance
    r rna abundance
    pmod proteinmodification
    var/variant
    frag/fragment
    loc/location
    bo/biologicalProcess
    pathology/path
    activity/act
    molecularActivity/ma
    translocation/tloc
    cellSecretion/sec
    cellSurfaceExpression/surf
    degradation/deg
    reacton/rxn
    fusion/fus
    
    """


    sourceid = allNodes[edge[0]]
    targetid = allNodes[edge[1]]

    interaction = edge[2] if edge[2] != None else ""

    elementsStr += "{ data: { source: '"+str(sourceid)+"', target: '"+str(targetid)+"', interact: '"+str(interaction)+"'}},"

elementsStr += "]\n}"


HTMLtop = """
<!DOCTYPE>
<!-- This code is for demonstration purposes only.  You should not hotlink to Github, Rawgit, or files from the Cytoscape.js documentation in your production apps. -->
<html>
  <head>
    <title>cose demo</title>
    <meta name="viewport" content="width=device-width, user-scalable=no, initial-scale=1, maximum-scale=1">
    <!-- For loading external data files -->
    <script src="cytoscape.js"></script> <!-- http://www.trauen-sich.net/js/ -->

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
        'target-arrow-color': '#9dbaea',
        'label': 'data(interact)'
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

with open(base+"network.html", 'w') as outHTML:

    outHTML.write( HTMLtop + elementsStr + HTMLbottom );