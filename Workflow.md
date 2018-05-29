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

Load the taxonomic file

```{r, message=FALSE}
taxonomy.raw = read.table("data/Saanich_cruise72_mothur_OTU_taxonomy.taxonomy", sep="\t", header=TRUE, row.names=1)
taxonomy= taxonomy.raw %>% 
  select(-Size) %>% 
  separate(Taxonomy, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
```

Change the name of the fasta alignment file 

```{r, message=FALSE}
alignment=read.fasta("data/mothur_intermediate_files/Saanich.final.opti_mcc.unique_list.0.03.rep.fasta")
names(alignment)=rownames(taxonomy)
```

Remove OTUS not assigned to a domain

```{r, message=FALSE}
taxonomy=subset(taxonomy,!Domain=="unknown")
```

Domain constrains

```{r, message=FALSE}
taxonomy[["Bacteria"]][taxonomy$Domain=="Bacteria"]=1
taxonomy[["Bacteria"]][taxonomy$Domain=="Archaea"]=0
```

Phylum constrains

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

Convert to fasta file

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

Now we are going to use the alignment of the reads to construct a phylogenetic tree of the reads.

To do so we are going to use the FastTree software, that you should already have installed on your computer. 

On mac, you need to open a terminal and acces the software commands where it is located (in my computer here: "Desktop/Programmes_Unix/FastTree")

The FastTree function requires input on (1) the model used (here gtr + cat 20), (2) the evenutal topological constrains (hre the file "Phylogenetic_Constrains.fasta" we produces before in R), (3) the alignment (here "Saanich.final.opti_mcc.unique_list.0.03.rep_Names_Modified.fasta") and a name for the output tree (here we used "Saanish_FastTree"). You need to give the exact adress of the different files in the command line  (here my files are loacted in "Desktop/Recherche/En_cours/workshopMicrobiome/Contenu/MyFiles/"). For further help on the the funciton, you can type: 

```{r, message=FALSE}
Desktop/Programmes_Unix/FastTree -h
```

Now we run fastTree on our data:

```{r, message=FALSE}
Desktop/Programmes_Unix/FastTree -gtr -cat 20 -constraints Desktop/Recherche/En_cours/workshopMicrobiome/Contenu/MyFiles/Phylogenetic_Constrains.fasta -nt Desktop/Recherche/En_cours/workshopMicrobiome/Contenu/MyFiles/Saanich.final.opti_mcc.unique_list.0.03.rep_Names_Modified.fasta > Desktop/Recherche/En_cours/workshopMicrobiome/Contenu/MyFiles/Saanish_FastTree
```

For the sake of comparison, we will laso construct the tree WITHOUT topological constrains. In this case just remove the parameters "constrains" 

```{r, message=FALSE}
#Desktop/Programmes_Unix/FastTree -gtr -cat 20 -nt Desktop/Recherche/En_cours/workshopMicrobiome/Contenu/MyFiles/Saanich.final.opti_mcc.unique_list.0.03.rep_Names_Modified.fasta > Desktop/Recherche/En_cours/workshopMicrobiome/Contenu/MyFiles/Saanish_FastTree_withoutConstrains
```

### 2.3. Tree Visualisation <a name="Tree-Visualisation"></a>

We now are coming back to R. We first need to load the trees we just constructed 

```{r, message=FALSE}
Tree=read.tree('My_outputs/Saanish_FastTree')
TreeNoC=read.tree('My_outputs/Saanish_FastTree_withoutConstrains')
```

We can now plot the tree, with colours for Domains. First creating a palette of colours for the Domains and assign each OTU a color based on that. 

```{r, message=FALSE}
Domains=unique(taxonomy$Domain)
paletteDomains=rainbow(length(Domains));names(paletteDomains)=Domains 
coloursDomains=paletteDomains[as.character(taxonomy$Domain)];names(coloursDomains)=rownames(taxonomy)
```

Then plot and write the tree in pdf format. 

```{r, message=FALSE}
pdf("My_outputs/Phylogenetic_tree_colouredby_Domains.pdf",width=15,height=15)
plot(Tree,type="fan",cex=.3,tip.color=coloursDomains[Tree$tip.label])
legend(1, 0, legend=names(paletteDomains),fill=paletteDomains, cex=2)
dev.off()
```

We can do the same for phylum instead of domain. 

```{r, message=FALSE}
Phylums=unique(taxonomy$Phylum)
palettePhylums=rainbow(length(Phylums));names(palettePhylums)=Phylums #define colours for taxonomic groups (Phylums here)
coloursPhylums=palettePhylums[as.character(taxonomy$Phylum)];names(coloursPhylums)=rownames(taxonomy) #assign colors to each OTU depending on its taxonomic group 
pdf("My_outputs/Phylogenetic_tree_colouredby_Phylum.pdf",width=15,height=15)
plot(Tree,type="fan",cex=.3,tip.color=coloursPhylums[Tree$tip.label])
legend(1, 1, legend=names(palettePhylums),fill=palettePhylums, cex=1)
dev.off()
```

It would be now interesting to check the impact of the topological contrains on the trees. To do so, we simply plot the tree obtained without constrains and observe the distribution of phylum in this tree. 

```{r, message=FALSE}
pdf("My_outputs/Phylogenetic_tree_colouredby_Phylum_NoConstrains.pdf",width=15,height=15)
plot(TreeNoC,type="fan",cex=.3,tip.color=coloursDomains[TreeNoC$tip.label])
legend(1, 1, legend=names(paletteDomains),fill=paletteDomains, cex=1)
dev.off()
```

Rapidly compare the two trees, what do you observe? 

### 2.3. Re-root the tree (Archaea Vs Bacteria) <a name="Tree-Rerooting"></a>

A phylogenetic tree is characterized by a root, which indicate the direction of evolution (from the root to the tip). For "prokariotes", we will assume to be on the branch separating bacteria from archaea. 

To do so, we first need to find this branch

```{r, message=FALSE}
Archaea=row.names(taxonomy)[taxonomy$Domain=="Archaea"]
MRCAnode=getMRCA(phy = Tree,tip = Archaea)
```

Then, we re-root the tre on this branch

```{r, message=FALSE}
TreeNewRoot=root(phy=Tree, node=MRCAnode,resolve.root = T)
```

and write the corresponding tree

```{r, message=FALSE}
write.tree(TreeNewRoot,'My_outputs/Saanish_FastTreeRooted')
```

We finally plot this tree and save the image. 

```{r, message=FALSE}
pdf("My_outputs/Phylogenetic_tree_colouredby_Domains_NewRoot.pdf",width=15,height=15)
plot(TreeNewRoot,type="fan",cex=.3,tip.color=coloursDomains[Tree$tip.label])
legend(1, 0, legend=names(paletteDomains),fill=paletteDomains, cex=2)
dev.off()
```

We now have a better tree. And we are going to use it to analyse microbiome composition patterns!

## 3. Diversity analysis in R <a name="Diversity-analysis-in-R"></a>

### 3.1. Basics <a name="Basics"></a>

Importing the data: OTU table, metadata and phylogenetic tree

The tree: 
```{r, message=FALSE}
Tree=read.tree('My_outputs/Saanish_FastTreeRooted')
```

The OTU table:
```{r, message=FALSE}
OTU = read.table("data/Saanich_cruise72_mothur_OTU_table.shared", sep="\t", header=TRUE, row.names=2)
OTU.clean = OTU %>%
  select(-label, -numOtus)
```

The taxonomic file:
```{r, message=FALSE}
taxonomy = read.table("data/Saanich_cruise72_mothur_OTU_taxonomy.taxonomy", sep="\t", header=TRUE, row.names=1)
taxonomy.clean = taxonomy %>% 
  select(-Size) %>% 
  separate(Taxonomy, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
```

The metadata:
```{r, message=FALSE}
metadata = read.table("data/Saanich_cruise72_metadata.txt", sep="\t", header=TRUE, row.names=1)
metadata=metadata[,c("Depth_m","PO4_uM","SiO2_uM","NO3_uM","NH4_uM","CH4_nM" ,"Salinity_PSU")]
```

We then Construct the phyloseq object:

First, we define parts of the phyloseq object.

```{r, message=FALSE}
OTU.clean.physeq = otu_table(as.matrix(OTU.clean), taxa_are_rows=FALSE)
tax.clean.physeq = tax_table(as.matrix(taxonomy.clean))
metadata.physeq = sample_data(metadata)
phylogeny.physeq=phy_tree(Tree)
```

and then assemble them

```{r, message=FALSE}
mothur = phyloseq(OTU.clean.physeq, tax.clean.physeq, metadata.physeq,phylogeny.physeq) 
mothur
```


### 3.2. Classic beta-diversity analysis <a name="Classics"></a>

We first compute classic metrics of beta-diversity

```{r, message=FALSE}
BC=vegdist(otu_table(mothur),method = "bray")
Jaccard=vegdist(otu_table(mothur),method = "jac")
UniFracBeta=UniFrac(mothur)
UniFracWBeta=UniFrac(mothur,weighted = T)
```

We can then plot them to see how much correlated they are
```{r, message=FALSE}
plot(Jaccard,UniFracBeta)
```

And then visualize the pattern of microbiome composition in relation to the metadata, for Unifrac:

```{r, message=FALSE}
ordi = ordinate(mothur, "PCoA", "unifrac", weighted=F)
plot_ordination(mothur, ordi, color="Depth_m")
```

or Bray Curtis: 

```{r, message=FALSE}
ordi = ordinate(mothur, "PCoA", "bray", weighted=F)
plot_ordination(mothur, ordi, color="Depth_m")
```

Importantly, we should also test the statistical relationship between microbiome beta-diversity and the metadata, using for example a PERMANOVA test: 

```{r, message=FALSE}
adonis(BC~Depth_m,data=data.frame(sample_data(mothur)))
adonis(Jaccard~Depth_m,data=data.frame(sample_data(mothur)))
adonis(UniFracBeta~Depth_m,data=data.frame(sample_data(mothur)))
```

	 
### 3.3. Phylogenetic beta-diversity <a name="Phylogenetic-Beta-diversity"></a>
	
### 3.4. Screening the phylogeentic scale: the BDTT analysis <a name="BDTT"></a>

### 3.5. Screening all the nodes of the phylogeny: PhyloFactor and Balance Trees <a name="creening-all-the-nodes"></a>

