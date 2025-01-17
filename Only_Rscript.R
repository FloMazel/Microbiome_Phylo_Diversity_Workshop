#' # Exploring the phylogenetic composition of microbiomes
#' 
#' Some of this code as been adapted from the [ECOSCOPE GitHub account](https://github.com/EDUCE-UBC) 
#' 
#' ## Table of contents
#' 
#' #### 1. [Prior to the workshop](#Prior)
#' 
#' #### 2. [Getting started with R ](#Introduction)
#' 
#' #### 3. [Building a microbial phylogenetic tree](#Building-a-microbial-phylogenetic-tree)
#' 
#' #### 4. [Diversity analysis in R ](#Diversity-analysis-in-R)
#' 
#' #### 5. [PhyloFactor](#screening-all-the-nodes)
#' 
#' ## 1. Prior to the workshop <a name="Prior"></a>
#' 
#' ### Setup 
#' 
#' Please come to the workshop with your laptop set up with the required
#' software and data files as described in our [setup instructions](https://github.com/FloMazel/Microbiome_Phylo_Diversity_Workshop/blob/master/SetUp.md).
#' 
#' ### Background 
#' 
#' Please read [Hallam SJ *et al*. 2017. Sci Data 4: 170158](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5663219/) "Monitoring microbial responses to ocean deoxygenation in a model oxygen minimum zone" to learn more about the data used in this workshop. You can also check out this [short video](https://drive.google.com/file/d/1s6V263-Vj2KQ6rExOBaYn-UDOERyQVGZ/view?usp=sharing) showing how the sampling was done!
#' 
#' 
#' ### Data description
#' The data in this workshop were collected as part of an on-going oceanographic time series program in Saanich Inlet, a seasonally anoxic fjord on the East coast of Vancouver Island, British Columbia 
#' 
#' 
#' ## 2. Getting started with R  <a name="Introduction"></a>
#' 
#' ### Installing and loading packages
#' 
#' At the beginning of every R script, you should have a dedicated space for loading R packages. R packages allow any R user to code reproducible functions and share them with the R community. Packages exist for anything ranging from microbial ecology to complex graphics to multivariate modeling and beyond. 
#' 
#' In this workshop, we will use several packages:
#' 
#' * tidyverse 
#' * ape
#' * seqinr
#' * vegan
#' * phyloseq
#' * betapart
#' * abind
#' * Matrix
#' * cowplot
#' * phylofactor
#' 
#' Note that `PhyloFactor` is not yet
#' available as a complete installable package on CRAN, so you will need to install
#' from source. To do this, execute the following commands in RStudio after you 
#' have installed all of the above packages.
#' 
## ------------------------------------------------------------------------
source("https://bioconductor.org/biocLite.R")
biocLite("ggtree")
biocLite("Biostrings")
install.packages('devtools')
devtools::install_github('reptalex/phylofactor', ref = '53ba93f02aee0e63c1e5a0ab2df0f5f976f08643')

#' 
#' 
#' Please MAKE SURE they are installed in your computer. If not, refer to [setup instructions](https://github.com/FloMazel/Microbiome_Phylo_Diversity_Workshop/blob/master/SetUp.Rmd)
#' 
#' ### Load the packages
#' 
## ---- message=FALSE------------------------------------------------------
library(seqinr)
library(ape)
library(vegan)
library(tidyverse)
library(phyloseq)
library(betapart)
library(abind)
library(tidyr)
library(Matrix)
library(cowplot)
library(phylofactor)

#' 
#' ### Setting your working directory 
#' 
#' Copy the workshop folder "Working directory" from GitHub on your computer
#' 
#' Then tell R that this is going to be the folder where we are going to work 
#' 
#' In Flo case: 
## ---- message=FALSE------------------------------------------------------
setwd("/Users/fmazel/Documents/GitHub/Microbiome_Phylo_Diversity_Workshop")

#' 
#' ### Importing custom R functions
#' 
#' We have written some custom R functions for use in this workshop. It is good
#' practice to define these functions in a separate file and then to import them
#' into your scripts, so that you don't duplicate the code. (Duplication, you 
#' may recall, can end up in mutation! ;) )
#' 
## ---- message=FALSE------------------------------------------------------
source("./R functions/BDTT_functions.R")

#' 
#' ## 3. Building a microbial phylogenetic tree <a name="Building-a-microbial-phylogenetic-tree"></a>
#' 
#' ### 3.1. Alignment <a name="Alignment"></a>
#' 
#' We are not going to cover this stage here, as it will be carried out in other workshops. 
#' Many software offer to align sequences, some of them are wrapped within [Mothur](https://www.mothur.org) or [Qiime2](https://qiime2.org)
#' For example you can look at this [tutorial](https://github.com/FloMazel/Microbiome_Phylo_Diversity_Workshop/blob/master/data/mothur_pipeline.html) producing the data we are going to use here. 
#' 
#' ### 3.2. Tree Building <a name="Tree-Building"></a>
#' 
#' Typical read length in microbiome studies are relatively short and have thus contains limited information to reconstruct phylogeentic trees, especially to reconstruct deep branches. To avoid biased phylogenies, we thus constrains deep branches to follow taxonomic classification as it is admitted that large taxonomic clades are monophyletic (informations based on longer sequences). The choice of the constrains is not easy. Here we will constain tree reconstruction by Domain (Bacteria/Archea) and phylums.
#' 
#' We will use FastTree to reconstruct the phylogeneic hypotheses. While fastTree is not the best software to reconstruct phylogenies (because it makes a lot of approximations), it has the huge avantage to be fast, which is often critical in microbial species, as there are a very high number of "species". 
#' 
#' #### Topological constrains
#' 
#' Load the taxonomic file
#' 
## ---- message=FALSE------------------------------------------------------
taxonomy.raw = read.table("data/Saanich_cruise72_mothur_OTU_taxonomy.taxonomy", sep="\t", header=TRUE, row.names=1)
taxonomy= taxonomy.raw %>% 
  select(-Size) %>% 
  separate(Taxonomy, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")

#' 
#' Change the name of the fasta alignment file 
#' 
## ---- message=FALSE------------------------------------------------------
alignment=read.fasta("data/mothur_intermediate_files/Saanich.final.opti_mcc.unique_list.0.03.rep.fasta")
names(alignment)=rownames(taxonomy)

#' 
#' Remove OTUS not assigned to a domain
#' 
## ---- message=FALSE------------------------------------------------------
taxonomy=subset(taxonomy,!Domain=="unknown")

#' 
#' Domain constrains
#' 
## ---- message=FALSE------------------------------------------------------
taxonomy[["Bacteria"]][taxonomy$Domain=="Bacteria"]=1
taxonomy[["Bacteria"]][taxonomy$Domain=="Archaea"]=0

#' 
#' Phylum constrains
#' 
## ---- message=FALSE------------------------------------------------------
Phylum=unique(taxonomy$Phylum)
Phylum=subset(Phylum,!Phylum%in%c("unknown_unclassified","Bacteria_unclassified", "Archaea_unclassified")) #remove this factor

for (i in Phylum)
  {
  taxonomy[[as.character(i)]][taxonomy$Phylum==i]=1
  taxonomy[[as.character(i)]][!taxonomy$Phylum==i]=0
  }

Constrains=taxonomy[,c("Bacteria",as.character(Phylum))] #keep only the constrains 
head(Constrains)

#' 
#' Convert to fasta file
#' 
## ---- message=FALSE------------------------------------------------------
sequences=list()
for (i in 1:dim(Constrains)[1]){sequences[[i]]=Constrains[i,]}
write.fasta(sequences, names=rownames(Constrains), file.out="My_outputs/Phylogenetic_Constrains.fasta", open = "w", nbchar = 60, as.string = FALSE)

#' 
#' Prune the alignment to sequences with assigned domains
#' 
## ---- message=FALSE------------------------------------------------------
alignment=alignment[rownames(Constrains)]
write.fasta(alignment, names=names(alignment), file.out="My_outputs/Saanich.final.opti_mcc.unique_list.0.03.rep_Names_Modified.fasta", open = "w", nbchar = 60, as.string = FALSE)

#' 
#' 
#' #### FastTree
#' 
#' Now we are going to use the alignment of the reads to construct a phylogenetic tree of the reads.
#' 
#' To do so we are going to use the FastTree software, that you should already have installed on your computer. 
#' 
#' On mac, you need to open a terminal and acces the software commands where it is located (in my computer here: "Desktop/Programmes_Unix/FastTree")
#' 
#' The FastTree function requires input on (1) the model used (here gtr + cat 20), (2) the evenutal topological constrains (hre the file "Phylogenetic_Constrains.fasta" we produces before in R), (3) the alignment (here "Saanich.final.opti_mcc.unique_list.0.03.rep_Names_Modified.fasta") and a name for the output tree (here we used "Saanish_FastTree"). You need to give the exact adress of the different files in the command line  (here my files are loacted in "Desktop/Recherche/En_cours/workshopMicrobiome/Contenu/MyFiles/"). For further help on the the funciton, you can type: 
#' 
## ---- message=FALSE------------------------------------------------------
#Desktop/Programmes_Unix/FastTree -h

#' 
#' Now we run fastTree on our data:
#' 
## ---- message=FALSE------------------------------------------------------
#Desktop/Programmes_Unix/FastTree -gtr -cat 20 -constraints Desktop/Recherche/En_cours/workshopMicrobiome/Contenu/MyFiles/Phylogenetic_Constrains.fasta -nt Desktop/Recherche/En_cours/workshopMicrobiome/Contenu/MyFiles/Saanich.final.opti_mcc.unique_list.0.03.rep_Names_Modified.fasta > Desktop/Recherche/En_cours/workshopMicrobiome/Contenu/MyFiles/Saanish_FastTree

#' 
#' For the sake of comparison, we will laso construct the tree WITHOUT topological constrains. In this case just remove the parameters "constrains" 
#' 
## ---- message=FALSE------------------------------------------------------
#Desktop/Programmes_Unix/FastTree -gtr -cat 20 -nt Desktop/Recherche/En_cours/workshopMicrobiome/Contenu/MyFiles/Saanich.final.opti_mcc.unique_list.0.03.rep_Names_Modified.fasta > Desktop/Recherche/En_cours/workshopMicrobiome/Contenu/MyFiles/Saanish_FastTree_withoutConstrains

#' 
#' ### 3.3. Tree Visualisation <a name="Tree-Visualisation"></a>
#' 
#' We now are coming back to R. We first need to load the trees we just constructed 
#' 
## ---- message=FALSE------------------------------------------------------
Tree=read.tree('My_outputs/Saanish_FastTree')
TreeNoC=read.tree('My_outputs/Saanish_FastTree_withoutConstrains')
Tree

#' 
#' We can now plot the tree, with colours for Domains. First creating a palette of colours for the Domains and assign each OTU a color based on that. 
#' 
## ---- message=FALSE------------------------------------------------------
Domains=unique(taxonomy$Domain)
paletteDomains=rainbow(length(Domains));names(paletteDomains)=Domains 
coloursDomains=paletteDomains[as.character(taxonomy$Domain)];names(coloursDomains)=rownames(taxonomy)

#' 
#' Then plot and write the tree in pdf format. 
#' 
## ---- message=FALSE------------------------------------------------------
pdf("My_outputs/Phylogenetic_tree_colouredby_Domains.pdf",width=15,height=15)
plot(Tree,type="fan",cex=.3,tip.color=coloursDomains[Tree$tip.label])
legend(1, 0, legend=names(paletteDomains),fill=paletteDomains, cex=2)
dev.off()

#' 
#' We can do the same for phylum instead of domain. 
#' 
## ---- message=FALSE------------------------------------------------------
Phylums=unique(taxonomy$Phylum)
palettePhylums=rainbow(length(Phylums));names(palettePhylums)=Phylums #define colours for taxonomic groups (Phylums here)
coloursPhylums=palettePhylums[as.character(taxonomy$Phylum)];names(coloursPhylums)=rownames(taxonomy) #assign colors to each OTU depending on its taxonomic group 
pdf("My_outputs/Phylogenetic_tree_colouredby_Phylum_with_Constrains.pdf",width=15,height=15)
plot(Tree,type="fan",cex=.3,tip.color=coloursPhylums[Tree$tip.label])
legend(1, 1, legend=names(palettePhylums),fill=palettePhylums, cex=1)
dev.off()

#' 
#' It would be now interesting to check the impact of the topological contrains on the trees. To do so, we simply plot the tree obtained without constrains and observe the distribution of phylum in this tree. 
#' 
## ---- message=FALSE------------------------------------------------------
pdf("My_outputs/Phylogenetic_tree_colouredby_Phylum_NoConstrains.pdf",width=15,height=15)
plot(TreeNoC,type="fan",cex=.3,tip.color=coloursPhylums[TreeNoC$tip.label])
legend(1, 1, legend=names(palettePhylums),fill=palettePhylums, cex=1)
dev.off()

#' 
#' Rapidly compare the two trees, what do you observe? 
#' 
#' ### 3.3. Re-root the tree (Archaea Vs Bacteria) <a name="Tree-Rerooting"></a>
#' 
#' A phylogenetic tree is characterized by a root, which indicate the direction of evolution (from the root to the tip). For "prokariotes", we will assume to be on the branch separating bacteria from archaea. 
#' 
#' To do so, we first need to find this branch
#' 
## ---- message=FALSE------------------------------------------------------
Archaea=row.names(taxonomy)[taxonomy$Domain=="Archaea"]
MRCAnode=getMRCA(phy = Tree,tip = Archaea)

#' 
#' Then, we re-root the tree on this branch
#' 
## ---- message=FALSE------------------------------------------------------
TreeNewRoot=root(phy=Tree, node=MRCAnode,resolve.root = T)

#' 
#' and write the corresponding tree
#' 
## ---- message=FALSE------------------------------------------------------
write.tree(TreeNewRoot,'My_outputs/Saanish_FastTreeRooted')

#' 
#' We finally plot this tree and save the image. 
#' 
## ---- message=FALSE------------------------------------------------------
pdf("My_outputs/Phylogenetic_tree_colouredby_Domains_NewRoot.pdf",width=15,height=15)
plot(TreeNewRoot,type="fan",cex=.3,tip.color=coloursDomains[Tree$tip.label])
legend(1, 0, legend=names(paletteDomains),fill=paletteDomains, cex=2)
dev.off()

#' 
#' We now have a better tree. And we are going to use it to document microbiome composition patterns!
#' 
#' ## 4. Diversity analysis in R <a name="Diversity-analysis-in-R"></a>
#' 
#' ### 4.1. Basics <a name="Basics"></a>
#' 
#' Importing the data: OTU table, metadata and phylogenetic tree
#' 
#' The tree: 
## ---- message=FALSE------------------------------------------------------
Tree=read.tree('My_outputs/Saanish_FastTreeRooted')

#' 
#' The OTU table:
## ---- message=FALSE------------------------------------------------------
OTU = read.table("data/Saanich_cruise72_mothur_OTU_table.shared", sep="\t", header=TRUE, row.names=2)
OTU.clean = OTU %>%
  select(-label, -numOtus)

#' 
#' The taxonomic file:
## ---- message=FALSE------------------------------------------------------
taxonomy = read.table("data/Saanich_cruise72_mothur_OTU_taxonomy.taxonomy", sep="\t", header=TRUE, row.names=1)
taxonomy.clean = taxonomy %>% 
  select(-Size) %>% 
  separate(Taxonomy, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")

#' 
#' The metadata:
## ---- message=FALSE------------------------------------------------------
metadata = read.table("data/Saanich_cruise72_metadata.txt", sep="\t", header=TRUE, row.names=1)
metadata=metadata[,c("Depth_m","PO4_uM","SiO2_uM","NO3_uM","NH4_uM","CH4_nM" ,"Salinity_PSU")]
head(metadata)

#' 
#' We then Construct the phyloseq object:
#' 
#' First, we define parts of the phyloseq object.
#' 
## ---- message=FALSE------------------------------------------------------
OTU.clean.physeq = otu_table(as.matrix(OTU.clean), taxa_are_rows=FALSE)
tax.clean.physeq = tax_table(as.matrix(taxonomy.clean))
metadata.physeq = sample_data(metadata)
phylogeny.physeq=phy_tree(Tree)

#' 
#' and then assemble them
#' 
## ---- message=FALSE------------------------------------------------------
saanish = phyloseq(OTU.clean.physeq, tax.clean.physeq, metadata.physeq,phylogeny.physeq) 
saanish

#' 
#' 
#' ### 4.2. Classic beta-diversity analysis <a name="Classics"></a>
#' 
#' We first compute classic metrics of beta-diversity
#' 
## ---- message=FALSE------------------------------------------------------
BC=vegdist(otu_table(saanish),method = "bray")
Jaccard=vegdist(otu_table(saanish),method = "jac")
UniFracBeta=UniFrac(saanish)

#' 
#' We can then plot them to see how much correlated they are
## ---- message=FALSE------------------------------------------------------
plot(Jaccard,UniFracBeta)

#' 
#' And then visualize the pattern of microbiome composition in relation to the metadata, for Unifrac:
#' 
## ---- message=FALSE------------------------------------------------------
ordi = ordinate(saanish, "PCoA", distance=UniFracBeta)
plot_ordination(saanish, ordi, color="Depth_m")

#' 
#' or Bray Curtis: 
#' 
## ---- message=FALSE------------------------------------------------------
ordi = ordinate(saanish, "PCoA", distance=BC)
plot_ordination(saanish, ordi, color="Depth_m")

#' 
#' Importantly, we should also test the statistical relationship between microbiome beta-diversity and the metadata, using for example a PERMANOVA test: 
#' 
## ---- message=FALSE------------------------------------------------------
adonis(BC~Depth_m,data=data.frame(sample_data(saanish)))
adonis(Jaccard~Depth_m,data=data.frame(sample_data(saanish)))
adonis(UniFracBeta~Depth_m,data=data.frame(sample_data(saanish)))
adonis(UniFracBeta~Depth_m+NO3_uM,data=data.frame(sample_data(saanish)))

#' 
#' Compare the correlation between environement and taxonomic or phylogenetic compositions
#' 
## ---- message=FALSE------------------------------------------------------
adonis(Jaccard~NO3_uM,data=data.frame(sample_data(saanish)))
adonis(UniFracBeta~NO3_uM,data=data.frame(sample_data(saanish)))

#' 
#' 
#' ### 4.3. Screening the phylogenetic scale <a name="BDTT"></a>
#' 
#' #### 4.3.1 Computing Beta-diversity
#' 
#' 
#' We now go a step ahead and screen the phylogentic scale to study the compositipon of the microbiomes. First, we check what are the values of this scale, i.e. the scale of divergence times, in substitution/site. 
#' 
## ---- message=FALSE------------------------------------------------------
Hnodes=getHnodes(Tree) 
hist(Hnodes,n=150,xlim=c(0,.5))

#' 
#' We can then define the set of discrete threshold used to build different set of microbial OTUs
#' 
## ---- message=FALSE------------------------------------------------------
slices=c(seq(from=0,to=0.3,by=0.025)) 

#' 
#' The idea the is to compute beta-doversity matrices for all these slices: 
#' 
#' We first extract a OTU table matrix from the phyloseq object
## ---- message=FALSE------------------------------------------------------
mat=t(as(otu_table(saanish), "matrix"))

#' 
#' and then run the analysis (for Jaccard and Bray Curtis) and explore the structure of the output (an "array") 
#' 
## ---- message=FALSE------------------------------------------------------
MultipleBetaJac=BDTT(similarity_slices=slices,tree=Tree,sampleOTUs=mat,onlyBeta=T,metric="jac")
class(MultipleBetaJac)
dim(MultipleBetaJac)
MultipleBetaJac
saveRDS(MultipleBetaJac,"My_outputs/Multiple_Resolution_Beta_Jaccard.RDS")  

#' 
#' 
#' Because it is slow to run (sorry for my code), we will not run bray curtis but directly load it from the BackUp
## ---- message=FALSE------------------------------------------------------
MultipleBetaBC=readRDS("My_outputs_BackUp/Multiple_Resolution_Beta_BrayCurtis.RDS")

#' 
#' #### 4.3.2. Statistical link to metadata
#' 
#' As for the simple beta-diversity tests, we then link the microbiome composition to our metadata. In our case, we need to construct models for each threshold independently:
#' 
#' We first prepare the result table: 
#' 
## ---- message=FALSE------------------------------------------------------
predictors=names(sample_data(saanish))
StatsRes=expand.grid(similarity_slices=as.character(slices),predictors=predictors,metric=c("Jac","BC"))
StatsRes[["F.Model"]]=StatsRes[["R2"]]=StatsRes[["Pr(>F)"]]=NA
head(StatsRes)

#' 
#' and then "fill" it without the models that we construct in loop: 
#' 
## ---- message=FALSE------------------------------------------------------
for (i in as.character(slices))
{
  for (j in predictors) 
  {
   res=unlist(adonis(formula = MultipleBetaJac[i,,]~data.frame(sample_data(saanish))[,j])$aov.tab[1,c(4,5,6)])
   StatsRes[(StatsRes$metric=="Jac")&(StatsRes$predictors==j)&(StatsRes$similarity_slices==i),4:6]=res
   res=unlist(adonis(formula = MultipleBetaBC[i,,]~data.frame(sample_data(saanish))[,j])$aov.tab[1,c(4,5,6)])
   StatsRes[(StatsRes$metric=="BC")&(StatsRes$predictors==j)&(StatsRes$similarity_slices==i),4:6]=res
   }
}

#' 
#' We can then plot the profiles of R2 along the phylogenetic time scale:
#' 
## ---- message=FALSE------------------------------------------------------
ggplot(aes(y=R2,x=similarity_slices,colour=predictors,group=factor(predictors)),data=StatsRes)+geom_point()+geom_line()+facet_wrap(~metric)
ggsave("My_outputs/BDTT_Jaccard_BC.pdf",height = 7,width = 13)

#' 
#' ## 5. Screening all the nodes of the phylogeny: PhyloFactor <a name="screening-all-the-nodes"></a>
#' 
#' ### 5.1. Loading data
#' 
#' We will load data the same data as above. PhyloFactor has some specific requirements for data format that will require further modifications, which are explained below.
#' 
#' The tree: 
## ---- message=FALSE------------------------------------------------------
Tree=read.tree('My_outputs/Saanish_FastTreeRooted')

#' 
#' The OTU table:
## ---- message=FALSE------------------------------------------------------
OTU = read.table("data/Saanich_cruise72_mothur_OTU_table.shared", sep="\t", header=TRUE, row.names=2)
OTU.clean = OTU %>%
  select(-label, -numOtus)

# phylofactor uses OTU tables in matrix format, with features as rows and samples as columns

OTU.matrix <- as.matrix(t(OTU.clean))

# we will need to subset the OTU table to only those tips found in the tree
OTU.matrix <- OTU.matrix[row.names(OTU.matrix) %in% Tree$tip.label,]

#' 
#' The taxonomic file:
## ---- message=FALSE------------------------------------------------------
taxonomy = read.table("data/Saanich_cruise72_mothur_OTU_taxonomy.taxonomy", sep="\t", header=TRUE, row.names=1)
taxonomy.clean = taxonomy %>% 
  select(-Size) %>% 
  separate(Taxonomy, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")

taxonomy <- taxonomy %>% select(-Size)
taxonomy$OTU_ID = row.names(taxonomy)
taxonomy <- taxonomy[Tree$tip.label, c('OTU_ID','Taxonomy')]
taxonomy[,2] <- as.character(taxonomy[,2])

#' 
#' The metadata:
## ---- message=FALSE------------------------------------------------------
metadata = read.table("data/Saanich_cruise72_metadata.txt", sep="\t", header=TRUE, row.names=1)
metadata=metadata[,c("Depth_m","PO4_uM","SiO2_uM","NO3_uM","NH4_uM","CH4_nM" ,"Salinity_PSU")]

# similarly, we need to subset metadata to just those samples in our data

metadata <- metadata[colnames(OTU.matrix),]

#' 
#' 
#' ### 5.2. Detecting a simulated variable on the bacterial phylogeny
#' 
#' Fundamentally, PhyloFactorization is a method for seeing how characteristics change across a phylogeny, and identifying the edges of the phylogeny where they change most significantly. 
#' 
#' We can demonstrate this by simulating some random data on top of our previously calculated microbial phylogeny, and adjusting the values for a few of the groups. You can imagine this simulated data as a characteristic of the microbes -- say, genome size, or cell diameter. 
#' 
#' Here, we're creating a new variable called `CellDiam`, and for each tip on our tree `Tree`, drawing a value from the log normal distribution using the function `rlnorm`. Then, we'll re-draw values for different subsets of the tree using values that are increased or decreased by 4-fold. 
#' 
## ------------------------------------------------------------------------

set.seed(1)
CellDiam <- rlnorm(length(Tree$tip.label))

# let's make Arcobacter and Actinobacteria different from the rest

# Grab a Boolean vector corresponding to the Arcobacter and Actinobacteria tip labels of the Tree object
arcobacter <- Tree$tip.label %in% row.names(taxonomy.clean[taxonomy.clean$Genus == 'Arcobacter',])
actinobacteria <- (Tree$tip.label %in% row.names(taxonomy.clean[taxonomy.clean$Phylum == 'Actinobacteria',]))

# Change those tips to have new values for CellDiam drawn from a different distribution
CellDiam[arcobacter] <- rlnorm(sum(arcobacter))*4
CellDiam[actinobacteria] <- rlnorm(sum(actinobacteria))/4

# Take the log of CellDiam variable for analysis
logCellDiam <- log(CellDiam)

#' 
#' Now we can PhyloFactor the tree to discover the branches that lead to the tips with the changed values. The function `twoSampleFactor` uses two-sample tests to identify edges separating two parts of the tree that differ with respect to a data vector -- in this case, the simulated `CellDiam` values.
#' 
## ------------------------------------------------------------------------
pf_twoSample <- twoSampleFactor(logCellDiam, Tree, nfactors=5)

#' 
#' We can look at a summary of the results by calling the object:
#' 
## ------------------------------------------------------------------------
pf_twoSample

#' 
#' Note that it returned 5 factors, but we only changed two clades on the tree! What happened?
#' 
#' PhyloFactor works by sequentially partitioning the tree from the most to the least important groups with respect to your given criteria. Note that the first two factors do in fact correspond to large groups, roughly matching the clades we adjusted! The remaining top five factors are individual tips that randomly obtained extreme draws from the simulated distribution.
#' 
## ------------------------------------------------------------------------
print(sum(actinobacteria))
print(sum(arcobacter))

#' 
#' And although statistical significance is still not well-understood for PhyloFactorization, you can see from the vector of p-values associated with these factors that there is a steep drop-off after the first two factors:
#' 
## ------------------------------------------------------------------------
pf_twoSample$pvals

#' 
#' We can also visualize these phylogenetic factors on the phylogeny itself using the `pf.tree` function in PhyloFactor, which creates a `ggtree` plot highlighting each returned factor:
#' 
## ---- message=FALSE------------------------------------------------------
cellDiam.tree <- pf.tree(pf_twoSample)

# Note that the object created by pf.tree has two parts: $legend, which contains 
# the colors used for each factor, and $ggplot, which has the ggplot object

cellDiam.tree$ggplot

#' 
#' Factor colors are written out as a separate legend element in the returned object:
#' 
## ------------------------------------------------------------------------
cellDiam.tree$legend

#' 
#' ### 5.3. Testing for phylogenetic groups correlated with metadata: linear regression
#' 
#' This will find the phylogenetic factors that vary most predictably with depth. By using the built-in choice function 'F', we are asserting that we want only the factors with the largest F score to be retained -- up to 5 factors. The `stop.early` parameter tells PhyloFactor to stop searching for new factors once it drops below a certain significance threshold (see the documentation for more on the heuristic used and its justification).
#' 
#' 
## ---- message=FALSE------------------------------------------------------
OTU.matrix <- OTU.matrix[Tree$tip.label,]
pf_depth_F <- PhyloFactor(OTU.matrix, Tree, X = metadata, frmla = Data~Depth_m, nfactors = 5, ncores=2, stop.early = T, choice='F')

#' 
#' (if you have an error on this function, you can load the backup stored object by uncommenting and running the line below)
#' 
## ---- message=FALSE------------------------------------------------------
# readRDS('./My_outputs_BackUp/pf_depth_F.RDS')

#' 
#' Looking at the pf_depth_F PhyloFactor object, we can see that it has primarily founds tips and small clades that are strongly correlated with depth.
#' 
## ------------------------------------------------------------------------
pf_depth_F

#' 
#' We can summarize the taxonomy of the first factor using `pf.taxa`:
## ------------------------------------------------------------------------
pf.taxa(pf_depth_F, taxonomy, factor=2)$group1

#' 
#' We see that it is a member of the Rhodobacteraceae, a photoheterotrophic lineage of Alphaproteobacteria. It would certainly make sense for this bacterium to be more abundant at shallower depths! We can use the `$data` element of the PhyloFactor summary object to plot the ILR-transformed relative abundances for this factor against depth, and see that it is in fact nicely correlated to depth:
#' 
## ------------------------------------------------------------------------
s_pf_depth.F.summary.1 <- summary(pf_depth_F, factor=1)

ggplot(s_pf_depth.F.summary.1$data, aes(x = Depth_m, y = Data)) + 
geom_point() + 
geom_smooth(method='lm')

#' 
#' 
#' We can also check out the taxonomy of the second factor, which actually comprises 19 tips in a single monophyletic clade.
#' 
## ------------------------------------------------------------------------
pf.taxa(pf_depth_F, taxonomy, factor=2)$group1

#' 
#' This clade of tips also belongs to the Rhodobacteraceae, and has a similar distribution:
#' 
## ------------------------------------------------------------------------
s_pf_depth.F.summary.2 <- summary(pf_depth_F, factor=2)

ggplot(s_pf_depth.F.summary.2$data, aes(x = Depth_m, y = Data)) + 
geom_point() + 
geom_smooth(method='lm')

#' 
#' ### 5.4. Partitioning the tree by total variance explained
#' 
#' The previous results were ranked according to the factors that had the strongest correlation, as measured by effect size (F). We can also rank factors according to the total proportion of variance explained, by setting `choice` equal to `'var'`.
#' 
## ---- message=FALSE------------------------------------------------------
pf_depth_var <- PhyloFactor(OTU.matrix, Tree, X = metadata, frmla = Data~Depth_m, nfactors = 5, ncores=2, stop.early = T, choice='var')

#' 
#' (if you have an error on this function, you can load the backup stored object by uncommenting and running the line below)
#' 
## ---- message=FALSE------------------------------------------------------
# readRDS('./My_outputs_BackUp/pf_depth_var.RDS')

#' 
## ------------------------------------------------------------------------
pf_depth_var

#' 
#' Note that when optimizing for total variance explained, PhyloFactor will tend to find bigger clades, as they represent a larger overall proportion of the data. The F-scores are lower in this case than for the smaller factors, which indicate that some members of the clade are not following the same linear relationship with depth.
#' 
#' We can see this, for examle, in the second PhyloFactor for the variance-paritioned tree. This factor represents a clade of Alphaproteobacteria, including the Rhodobacteraceae we found above:
#' 
## ------------------------------------------------------------------------
pf.taxa(pf_depth_var, taxonomy, factor=2)$group1

#' 
#' As you can see by plotting the ILR abundances of this entire clade against depth, the correlation is not as tight across the full range. Other bacteria in this clade may be especially more prominant at the very shallowest depths:
#' 
## ------------------------------------------------------------------------
s_pf_depth.var.summary.2 <- summary(pf_depth_var, factor=2)

ggplot(s_pf_depth.var.summary.2$data, aes(x = Depth_m, y = Data)) + 
geom_point() + 
geom_smooth(method='lm')

#' 
#' 
#' ### 5.5. Comparing phylofactorization with different independent variables
#' 
#' Let's combine all we've practiced so far to visualize using PhyloFactorization how different clades of the bacterial phylogeny respond to changes in different variables in our dataset. 
#' 
#' We can use the `cowplot` package to combine different ggplot elements to simultaneously illustrate the ILR-values for individual PhyloFactors and their location on the phylogeny. For this, we'll just take the top three factors 
#' 
#' First, we'll write a short helper function to format the individual inset XY plots.
#' 
## ------------------------------------------------------------------------

inset_plot <- function(pf, factor, variable, color) {
  pf_s <- summary(pf, factor=factor)
  
  p <- ggplot(pf_s$data, aes_string(x = variable, y = 'Data')) + 
        geom_point(color = color) + 
        geom_abline(slope = coef(pf$models[[factor]])[2],
                    intercept = coef(pf$models[[factor]])[1])
  
  return(p)
}

#' 
#' Now we can perform a phylofactorization of the bacterial phylogeny with respect to different variables of interest. We'll compare the bacterial groups responding to nitrate (NO3) and phosphate (PO4), two important and frequently limiting nutrients in oligotrophic environments. Let's grab the top three factors for each, according to total variance explained.
#' 
## ---- message=FALSE------------------------------------------------------
pf_NO3_var <- PhyloFactor(OTU.matrix, Tree, X = metadata, frmla = Data~NO3_uM, nfactors = 3, ncores=2, choice='var')

pf_PO4_var <- PhyloFactor(OTU.matrix, Tree, X = metadata, frmla = Data~PO4_uM, nfactors = 3, ncores=2, choice='var')

#' 
#' 
#' (if you have an error on this function, you can load the backup stored object by uncommenting and running the lines below)
#' 
## ---- message=FALSE------------------------------------------------------
# readRDS('./My_outputs_BackUp/pf_NO3_var.RDS')
# readRDS('./My_outputs_BackUp/pf_PO4_var.RDS')

#' 
#' 
#' For each PhyloFactor object, we can use the `cowplot` library to combine separate `ggplot2` panels into a publication-ready image, showing both individual factors locations on the tree, and their response (in ILR-transformed space) to each predictor variable.
#' 
## ---- message=FALSE------------------------------------------------------
# Draw the PhyloFactors on the phylogeny
NO3_tree = pf.tree(pf_NO3_var)

# Make inset plots for each factor
g1_NO3 = inset_plot(pf_NO3_var,                  # The PhyloFactor object we just calculated
                    1,                           # The first factor
                    'NO3_uM',                    # Plotting against the NO3_uM variable
                    NO3_tree$legend$colors[1])   # Use the color for this factor from the tree
g2_NO3 = inset_plot(pf_NO3_var, 2, 'NO3_uM', NO3_tree$legend$colors[2])
g3_NO3 = inset_plot(pf_NO3_var, 3, 'NO3_uM', NO3_tree$legend$colors[3])

# Draw the inset plots composed with the PhyloFactor phylogeny
ggdraw() +
  draw_plot(NO3_tree$ggplot, 
            x = 0, y = 0.05, width = .7, height = .7,scale = 1.8) +
  draw_plot(g1_NO3, x = .7, y = .66, width = .3, height = .33) +
  draw_plot(g2_NO3, x = .7, y = .33, width = .3, height = .33) +
  draw_plot(g3_NO3, x = .7, y = 0, width = .3, height = .33) +
  draw_plot_label(label = c("Phylogeny", "A", "B","C"), size = 15,
                  x = c(0, 0.67, .67,.67), y = c(1, 1, 0.66,0.33))


#' 
## ---- message=FALSE------------------------------------------------------
# Repeat for the phosphate PhyloFactors
PO4_tree = pf.tree(pf_PO4_var)

g1_PO4 = inset_plot(pf_PO4_var, 1, 'PO4_uM', PO4_tree$legend$colors[1])
g2_PO4 = inset_plot(pf_PO4_var, 2, 'PO4_uM', PO4_tree$legend$colors[2])
g3_PO4 = inset_plot(pf_PO4_var, 3, 'PO4_uM', PO4_tree$legend$colors[3])

ggdraw() +
  draw_plot(PO4_tree$ggplot, 
            x = 0, y = 0.05, width = .7, height = .7,scale = 1.8) +
  draw_plot(g1_PO4, x = .7, y = .66, width = .3, height = .33) +
  draw_plot(g2_PO4, x = .7, y = .33, width = .3, height = .33) +
  draw_plot(g3_PO4, x = .7, y = 0, width = .3, height = .33) +
  draw_plot_label(label = c("Phylogeny", "A", "B","C"), size = 15,
                  x = c(0, 0.67, .67,.67), y = c(1, 1, 0.66,0.33))


#' 
#' ### 5.6. Other visualizations using PhyloFactor: phylogenetic heatmaps
#' 
#' Heatmaps can help give us a better sense of how the data are actually distributed across samples and the tree. We can use the `pf.heatmap` function to draw both the tree (and highlight the relevant PhyloFactor partitions) and the (ILR-transformed) raw data for each tip:
#' 
## ---- message=FALSE------------------------------------------------------
pf.heatmap(pf_PO4_var, 
           factors=1:3, 
           column.order = order(metadata$PO4_uM),
           low='purple',
           high='yellow')

#' 
#' Heatmaps can also be used to illustrate the values predicted across the tree, using the linear models fit by PhyloFactorization above:
#' 
## ---- message=FALSE------------------------------------------------------

preds <- predict(pf_PO4_var)

pf.heatmap(pf_PO4_var, 
           factors=1:3, 
           Data=preds,
           column.order = order(metadata$PO4_uM),
           low='purple',
           high='yellow')

#' 
