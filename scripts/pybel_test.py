import json
import requests
import time

import networkx as nx

import pybel
from pybel.constants import *
import pybel_tools
from pybel_tools.visualization import to_jupyter

base = "/mnt/c/ownCloud/data/bcn/"
inputFile = base + "all/CV-IPN-Endothelial cell activation1.0.sif"
inputFile = base + "manual/CV-IPN-Endothelial cell activation1.0.jgf"
file = open(inputFile, 'r')

res = json.load(file)

def get_citation(evidence):
    return {
        CITATION_NAME: evidence['citation']['name'],
        CITATION_TYPE: evidence['citation']['type'],
        CITATION_REFERENCE: evidence['citation']['id']
    }


annotation_map = {
    'tissue': 'Tissue',
    'disease': 'Disease',
    'species_common_name': 'Species'
}


species_map = {
    'human': '9606',
    'rat': '10116',
    'mouse': '10090'
}


annotation_value_map = {
    'Species': species_map
}

graph = pybel.BELGraph()
parser = pybel.parser.BelParser(graph)

for edge in res['graph']['edges']:
    for evidence in edge['metadata']['evidences']:
        if 'citation' not in evidence or not evidence['citation']:
            continue

        d = {}

        if 'biological_context' in evidence:
            annotations = evidence['biological_context']

            if annotations['tissue']:
                d['Tissue'] = annotations['tissue']

            if annotations['disease']:
                d['Disease'] = annotations['disease']

            if annotations['species_common_name']:
                d['Species'] = species_map[annotations['species_common_name'].lower()]

        print(d)
        print(evidence)

        #parser.control_parser.annotations.update(d)

    print(edge)
    print(edge['metadata']['evidences'])

    parser.control_parser.clear()
    parser.control_parser.citation = {
        CITATION_NAME: "MJ",
        CITATION_TYPE: "custom",
        CITATION_REFERENCE: "MJ"
    }
    #parser.control_parser.evidence = evidence['summary_text']

    bel = '{source} {relation} {target}'.format_map(edge)
    try:
        tokenized = parser.parseString(bel)
        print(tokenized)
    except Exception as e:

        print(e, bel)
        exit(0)