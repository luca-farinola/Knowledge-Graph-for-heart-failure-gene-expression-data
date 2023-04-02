# -*- coding: utf-8 -*-
"""
Created on Sun Feb 12 11:12:38 2023    
@author: lukfa
"""
import pandas as pd
import networkx as nx
import os
from rdflib import URIRef,Literal, Namespace


mynamespace = Namespace('http://lucafarinola.com/my_vocabulary#')
ns = Namespace("http://bio2rdf.org/")


#put your directory 
os.chdir('C:/Users/lukfa/OneDrive/Desktop/Diffexpr')

# Load the differential expression table
genes = pd.read_csv("data/MAGNET_GeneExpressionData_CPM_19112020.txt", sep = '\t')
diffexpr = pd.read_csv("files/dcmvscontroll_genelist.txt", sep = '\t',index_col=0)
genes = genes.loc[genes['EnsemblGeneID'].isin(diffexpr.index)]
genes.index = diffexpr['hgnc_id']
genes = genes.T
genes.drop(index = 'EnsemblGeneID',inplace = True)
genes.columns= genes.columns.str.lower()

# Calculate the pairwise Pearson correlation between the genes
correlation_matrix = genes.astype(float).corr(method='spearman')

data = correlation_matrix.to_numpy()

nodes = correlation_matrix.columns

# Create a NetworkX graph object
G = nx.Graph()
G.add_nodes_from(nodes)

# Set the threshold for retaining only the most highly correlated gene pairs
threshold = 0.7

for i, node1 in enumerate(nodes):
    for j, node2 in enumerate(nodes):
        if i < j and abs(data[i][j]) > threshold:
            G.add_edge(node1, node2, weight=data[i][j])
            
             
import rdflib

nodes_with_edges = [node for node in G.nodes() if len(G.edges(node)) > 0]

#add info with biotordf
biotordf = rdflib.Graph()
my = rdflib.Graph()

for i in nodes_with_edges: 
    if str(i) != 'nan':
        biotordf.parse('http://bio2rdf.org/' + str(i))
        
biotordf.serialize(destination='files/biotordf.ttl', format='turtle') 
   
#schema for biotordf
ns1 = 'http://bio2rdf.org/hgnc_vocabulary:'
ns2 = 'http://purl.org/dc/terms/'
ns3 = 'http://rdfs.org/ns/void#'
ns4 = 'http://bio2rdf.org/bio2rdf_vocabulary:'
rdfs = 'http://www.w3.org/2000/01/rdf-schema#'
xsd = 'http://www.w3.org/2001/XMLSchema#'
   
#shema for my KG
hgnc = Namespace('http://identifiers.org/hgnc:')
fam = Namespace('http://identifiers.org/hgnc.family/')
biolink = Namespace("https://w3id.org/biolink/vocab/")
pubmed = Namespace('https://pubmed.ncbi.nlm.nih.gov/')
omim = Namespace('https://www.omim.org/entry/')


# for omim or https://www.omim.org/entry/ ? 
# http://purl.obolibrary.org/obo/OMIM_
     
for s, p, o in biotordf:
    a = URIRef('http://www.w3.org/1999/02/22-rdf-syntax-ns#type')
    g = URIRef(str(hgnc + s[24:])) 
    #change predicates  
    if o == URIRef(ns1 + 'Gene-Symbol'):
        my.add((g,a,URIRef(biolink + 'Gene')))
    elif p == URIRef(ns1 + 'chromosome'): 
        pred = URIRef(biolink + 'expressed_in')
        my.add((g,pred,Literal(o)))
    elif p == URIRef(ns1 + 'approved-symbol'): 
        pred = URIRef(biolink + 'full_name')
        my.add((g,pred,Literal(o)))
    elif p == URIRef(ns1 + 'gene-family-description'): 
        pred = URIRef(biolink + 'subclass_of')
        my.add((g,pred,Literal(o)))
    elif p == URIRef(ns1 + 'x-pubmed'):
        pred = URIRef(biolink + 'is_assessed_by')
        obj = URIRef(pubmed + o[26:])
        my.add((g,pred,obj))
    elif p == URIRef(ns1 + 'x-omim'): 
        pred = URIRef(biolink + 'has_phenotype')
        obj = URIRef(omim + o[26:])
        my.add((g,pred,obj))
    elif p == URIRef(ns1 + 'xref'):
        if o[:26]  == 'http://cancer.sanger.ac.uk':
            pred = URIRef(biolink + 'risk_affected_by')
            my.add((g,pred,o))


# Define the edges from the gene coexpression network
for edge in G.edges(data=True):
    subject_uri = URIRef('http://identifiers.org/' + str(edge[0]))
    object_uri = URIRef('http://identifiers.org/' + str(edge[1]))
    correlates = "positively_correlated_with"
    anticorrelates = "negatively_correlated_with"
    if edge[2]['weight'] > 0:
        my.add((subject_uri, URIRef(biolink + correlates), object_uri))
    if edge[2]['weight'] < 0:
        my.add((subject_uri, URIRef(biolink + anticorrelates), object_uri))
        

len(my)
diffexpr.index = diffexpr['hgnc_id'].str.lower()
# add foldchange and pvalues
for index,row in diffexpr.iterrows():
    if str(index) in nodes_with_edges :
            subject_uri = URIRef('http://identifiers.org/'  + str(index)) 
            fc = diffexpr.at[index,'logFC']
            pval = diffexpr.at[index,'adj.P.Val']
            upfoldchange = "increases_amount_or_activity_of"
            downfoldchange = "decreases_amount_or_activity_of"
            haspval = "adjusted_p_value" # am I allowed ? it's biolink edge property ?
            my.add((subject_uri, URIRef(biolink + haspval),Literal(pval)))
            if fc > 0 :
                my.add((subject_uri, URIRef(biolink + upfoldchange),Literal(fc)))
            elif fc < 0 :
                my.add((subject_uri, URIRef(biolink + downfoldchange),Literal(fc)))

ontologies = pd.read_csv("files/geneontology_DCM.txt", sep = ',')
go = Namespace('http://amigo.geneontology.org/amigo/term/GO:')


for index,row in ontologies.iterrows():
    GO = str(ontologies.at[index,'ID'])
    HGNC = str(ontologies.at[index,'hgnc_id']).lower()
    ont_type = str(ontologies.at[index,'name_1006'])
    ont_class = str(ontologies.at[index,'namespace_1003'])
    pval = ontologies.at[index,'p.adjust']
    if HGNC in nodes_with_edges :
            gene = URIRef('http://identifiers.org/' + HGNC)  
            ontology = URIRef(go + GO[3:])
            a = URIRef('http://www.w3.org/1999/02/22-rdf-syntax-ns#type')
            if len(GO) > 3: # avoid nan 
                my.add((gene,URIRef(biolink + 'participates_in'),ontology ))
                my.add((ontology, a ,URIRef(biolink + 'OntologyClass')))
                my.add((ontology,URIRef(biolink + 'part_of'), Literal(ont_class)))
                my.add((ontology, URIRef(biolink + 'associated_with') , Literal(ont_type)))
                my.add((ontology, URIRef(biolink + 'adjusted_p_value') , Literal(pval)))
                                


#add stuff about GO 
my.bind("rdfs", rdfs)
my.bind("hgnc", hgnc)
my.bind("pubmed", pubmed)
my.bind("omim", omim)
my.bind("biolink", biolink)
my.bind("go", go)
my.serialize(destination='files/prova.ttl', format='turtle')

#networkX 
from rdflib.extras.external_graph_libs import rdflib_to_networkx_digraph
from rdflib.extras.external_graph_libs import rdflib_to_networkx_multidigraph
from operator import itemgetter

mdg = rdflib_to_networkx_multidigraph(my)
dg = rdflib_to_networkx_digraph(my)

mdg.number_of_nodes() 

degrees = [val for (node, val) in mdg.degree()]  # gets the degrees of each node in the network
sum_degrees = sum(degrees) # sums up all the degrees of all nodes in the network
avg_degree_g = sum_degrees / mdg.number_of_nodes() # calculate average degree of network by dividing sum of all node degrees and dividing by number of nodes
avg_degree_g 

degree_dict = dict(mdg.degree(mdg.nodes()))
nx.set_node_attributes(mdg, degree_dict, 'degree')

sorted_degree = sorted(degree_dict.items(), key=itemgetter(1), reverse=True)

print("Top 50 nodes by degree:")
for d in sorted_degree[:50]:
    print(d)

results = nx.degree_centrality(mdg)

for node_id, rank in sorted(results.items(), key=lambda item: item[1], reverse=True)[:10]:
        print("{:6.3f} {}".format(rank, node_id))

# Rule Mining 

query = """
SELECT ?subject ?predicate ?object WHERE { 
  ?subject ?predicate ?object 
}
"""

results = my.query(query)
# Write the results to a TSV file
with open("my_kg.tsv", "w") as f:
    for row in results:
        f.write("\t".join(str(x) for x in row) + "\n")

        
import os
os.system('java -jar C:/Users/lukfa/OneDrive/Desktop/KG/lab8/lab_rule_mining/amie_plus.jar my_kg.tsv > amie_result.tsv')


# KG exploration 

pubmed = []

for s, p, o in my:
    if p == URIRef(biolink + 'is_assessed_by'):
        o = str(o)
        pubmed.append(o)


def common(strings):
    string_counts = {}
    for string in strings:
        if string in string_counts:
            string_counts[string] += 1
        else:
            string_counts[string] = 1
    return string_counts

pubmed = common(pubmed) 

for key in pubmed: 
    if pubmed[key] > 1 : 
        print(key,pubmed[key])
        
