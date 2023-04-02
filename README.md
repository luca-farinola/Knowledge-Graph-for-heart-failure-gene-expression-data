  <h3 align="center">KG for heart failure gene expression data</h3>

  <p align="center">
    Project for Building and Mining Knowledge Graph 
  </p>
</div>

<!-- ABOUT THE PROJECT -->
## About The Project

In this project gene expression data are used to construct a Knowledge graph. 

Bird eye view of the project:
* Differential expression analysis with limma (R)
* Co-epression Networks wih NetworkX (Python)
* Knowledge graph construction with rdflib (Python)
* Queries to retrieve interesting genes/biological processes (SPARQL) 

Main sources of information are : Biomart, Bio2RDF,Coexpression matrix and Differential expression values.
Genes are connected to Molecular functions, Biological Processes, Cellular components through EnrichGO but also to pvalues, foldchange
or other URL such as COSMIC or PubMed articles. 

<img src="Images/KG_metagraph.png" widht = "400" height= "400">

<!-- DATA -->
### Data

Data directory is needed to run the files. Those data include gene expression count data, 
metadata with patient information, gene lenght for FPKM normalizzation. Those data are uploaded in zenodo at link below.  

</a>
<a href="https://doi.org/10.5281/zenodo.7790931">
        <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.7790931.svg" alt="DOI">
    </a>
