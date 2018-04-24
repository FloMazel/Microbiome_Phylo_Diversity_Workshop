# Exploring the phylogenetic composition of microbiomes


## Table of contents


#### 1. [Introduction](#Introduction)

##### 2. [Building a microbial phylogenetic tree](#Building-a-microbial-phylogenetic-tree)

##### 2.1. [Alignment](#Alignment)

##### 2.2. [Tree Building](#Tree-Building)

#### 3. [Diversity analysis in R ](#Diversity-analysis-in-R)

##### 3.1. [Basics ](#Basics)

##### 3.2. [Taxonomic Beta-diversity](#Taxonomic-Beta-diversity)

##### 3.3.[Phylogenetic Beta-diversity](#Phylogenetic-Beta-diversity)

##### 3.4.[Screening the phylogeentic scale: the BDTT analysis](#BDTT)

##### 3.5 [Screening all the nodes of the phylogeny: PhyloFactor and Balance Trees](#Screening-all-the-nodes)



## 1. Introduction <a name="Introduction"></a>

Alpha, Beta and Gamma diversity
Microbial units and the need of a phylogenetic framework


## 2. Building a microbial phylogenetic tree <a name="Building-a-microbial-phylogenetic-tree"></a>

### 2.1. Alignment <a name="Alignment"></a>

We are not going to cover this stage here, as it will be carried out in other workshops. 
It can be done within *Mothur or *Qiime2


### 2.2. Tree Building <a name="Tree-Building"></a>

Using FastTree (pyhton) and R here

Typical read length in microbiome studies are relatively short and have thus contains limited information to reconstruct phylogeentic trees, especially to reconstruct deep branches. To avoid biased phylogenies, we thus constrains deep branches to follow taxonomic classificaiton as we know that large taxonomic clades are made to be monphyletic (informations based on longer sequences). The choice of the constrains is not easy. Here we will constain tree reconstruction by Domain (Bacteria/Archea) and phylums.

We will use FastTree to reconsttuct the phylogeneic hypotheses. While fastTree is not the best software to reconstruct phylogenies (because it made a lot of approaximations), it has the huge avantage to be fast, which is often critical in microbial species, as there are a very high number of species. 

#### Topological constrains

Using R 

See corresponding code for now


#### FastTree

Using the terminal 
	

## 3. Diversity analysis in R <a name="Diversity-analysis-in-R"></a>

### 3.1. Basics <a name="Basics"></a>

Importing the data: OTU table, metadat and phylogenetic tree

### 3.2. Taxonomic Beta-diversity <a name="Taxonomic-Beta-diversity"></a>
	 
### 3.3. Phylogenetic beta-diversity <a name="Phylogenetic-Beta-diversity"></a>
	
### 3.4. Screening the phylogeentic scale: the BDTT analysis <a name="BDTT"></a>

### 3.5. Screening all the nodes of the phylogeny: PhyloFactor and Balance Trees <a name="creening-all-the-nodes"></a>



## Introduction <a name="Introduction"></a>

