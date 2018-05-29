# Exploring the phylogenetic composition of microbiomes

Some of this code as been adapted from the [ECOSCOPE GitHub account](https://github.com/EDUCE-UBC) 

## Table of contents

#### 1. [Prior to the workshop](#Prior)

#### 2. [Introduction](#Introduction)

#### 3. [Building a microbial phylogenetic tree](#Building-a-microbial-phylogenetic-tree)

#### 4. [Diversity analysis in R ](#Diversity-analysis-in-R)

## 1. Prior to the workshop <a name="Prior"></a>

### Setup 

Please come to the workshop with your laptop set up with the required
software and data files as described in our [setup instructions](https://github.com/FloMazel/Microbiome_Phylo_Diversity_Workshop/blob/master/SetUp.md).

### Background 

Please read [Hallam SJ *et al*. 2017. Sci Data 4: 170158](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5663219/) "Monitoring microbial responses to ocean deoxygenation in a model oxygen minimum zone" to learn more about the data used in this workshop. You can also check out this [short video](https://drive.google.com/file/d/1s6V263-Vj2KQ6rExOBaYn-UDOERyQVGZ/view?usp=sharing) showing how the sampling was done!


### Data description
The data in this workshop were collected as part of an on-going oceanographic time series program in Saanich Inlet, a seasonally anoxic fjord on the East coast of Vancouver Island, British Columbia 


## 2. Introduction <a name="Introduction"></a>

### 2.1. Theory
General concepts: Alpha, Beta and Gamma diversity
Microbial units and the need of a phylogenetic framework

### 2.2. Getting started with R 

#### Installing and loading packages

At the beginning of every R script, you should have a dedicated space for loading R packages. R packages allow any R user to code reproducible functions and share them with the R community. Packages exist for anything ranging from microbial ecology to complex graphics to multivariate modeling and beyond. 

In this workshop, we will use several packages 
XXX

```{r, message=FALSE}
if (!require("seqinr")) install.packages("seqinr")
if (!require("ape")) install.packages("ape")
if (!require("vegan")) install.packages("vegan")
```

#### Setting your working directory 

Copy the workshop folder "Working directory" from GitHub on your computer

Then tell R that this is going to be the folder where we are going to work 

In my case: 
```{r, message=FALSE}
setwd("/Users/fmazel/Desktop/Recherche/En_cours/workshopMicrobiome/WorkingDirectory")
```

## 3. Building a microbial phylogenetic tree <a name="Building-a-microbial-phylogenetic-tree"></a>

### 2.1. Alignment <a name="Alignment"></a>

We are not going to cover this stage here, as it will be carried out in other workshops. 
Many software offer to align sequences, some of them are wrapped within [Mothur](https://www.mothur.org) or [Qiime2](https://qiime2.org)
For example you can look at this [tutorial](https://github.com/FloMazel/Microbiome_Phylo_Diversity_Workshop/blob/master/data/mothur_pipeline.html) producing the data we are going to use here. 

### 2.2. Tree Building <a name="Tree-Building"></a>

Typical read length in microbiome studies are relatively short and have thus contains limited information to reconstruct phylogeentic trees, especially to reconstruct deep branches. To avoid biased phylogenies, we thus constrains deep branches to follow taxonomic classification as it is admitted that large taxonomic clades are monophyletic (informations based on longer sequences). The choice of the constrains is not easy. Here we will constain tree reconstruction by Domain (Bacteria/Archea) and phylums.

We will use FastTree to reconsttuct the phylogeneic hypotheses. While fastTree is not the best software to reconstruct phylogenies (because it makes a lot of approximations), it has the huge avantage to be fast, which is often critical in microbial species, as there are a very high number of species. 

#### Topological constrains

##### Load the taxonomic file

```{r, message=FALSE}
taxonomy.raw = read.table("data/Saanich_cruise72_mothur_OTU_taxonomy.taxonomy", sep="\t", header=TRUE, row.names=1)
taxonomy= taxonomy.raw %>% 
  select(-Size) %>% 
  separate(Taxonomy, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
```

##### Change the name of the fasta alignment file 

```{r, message=FALSE}
alignment=read.fasta("data/mothur_intermediate_files/Saanich.final.opti_mcc.unique_list.0.03.rep.fasta")
names(alignment)=rownames(taxonomy)
```

##### Remove OTUS not assigned to a domain

```{r, message=FALSE}
taxonomy=subset(taxonomy,!Domain=="unknown")
```

##### Domain constrains

```{r, message=FALSE}
taxonomy[["Bacteria"]][taxonomy$Domain=="Bacteria"]=1
taxonomy[["Bacteria"]][taxonomy$Domain=="Archaea"]=0
```

##### Phylum constrains

```{r, message=FALSE}
Phylum=unique(taxonomy$Phylum)
Phylum=subset(Phylum,!Phylum=="unknown_unclassified") #remove this factor

for (i in Phylum)
  {
  taxonomy[[as.character(i)]][taxonomy$Phylum==i]=1
  taxonomy[[as.character(i)]][!taxonomy$Phylum==i]=0
  }

Constrains=taxonomy[,c("Bacteria",as.character(Phylum))] #keep only the constrains 
```

##### Convert to fasta file

```{r, message=FALSE}
sequences=list()
for (i in 1:dim(Constrains)[1]){sequences[[i]]=Constrains[i,]}
write.fasta(sequences, names=rownames(Constrains), file.out="My_outputs/Phylogenetic_Constrains.fasta", open = "w", nbchar = 60, as.string = FALSE)
```

##### Prune the alignment to sequences with assigned domains

```{r, message=FALSE}
alignment=alignment[rownames(Constrains)]
write.fasta(alignment, names=names(alignment), file.out="My_outputs/Saanich.final.opti_mcc.unique_list.0.03.rep_Names_Modified.fasta", open = "w", nbchar = 60, as.string = FALSE)
```


#### FastTree

Using the terminal 

See corresponding code (for now)

## 3. Diversity analysis in R <a name="Diversity-analysis-in-R"></a>

### 3.1. Basics <a name="Basics"></a>

Importing the data: OTU table, metadata and phylogenetic tree

### 3.2. Taxonomic Beta-diversity <a name="Taxonomic-Beta-diversity"></a>
	 
### 3.3. Phylogenetic beta-diversity <a name="Phylogenetic-Beta-diversity"></a>
	
### 3.4. Screening the phylogeentic scale: the BDTT analysis <a name="BDTT"></a>

### 3.5. Screening all the nodes of the phylogeny: PhyloFactor and Balance Trees <a name="creening-all-the-nodes"></a>

