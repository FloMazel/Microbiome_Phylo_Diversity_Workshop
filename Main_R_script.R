################################################################################
# load and install packages
################################################################################

# if the named package is not installed then install it
if (!require("seqinr")) install.packages("seqinr")
if (!require("ape")) install.packages("ape")
if (!require("vegan")) install.packages("vegan")

#load packages
library(seqinr)
library(ape)

#Put here your working directory
setwd("/Users/fmazel/Desktop/Recherche/En_cours/workshopMicrobiome/WorkingDirectory")

################################################################################
# Part I: Construct the phylogenetic Tree
################################################################################

# Part 1.2. Construct the constrains file based on taxonomic information 
################################################################################

# load the taxonomic file
taxonomy=read.table("data/Saanich_cruise72_mothur_OTU_taxonomy.modified.txt",row.names=1,header=T)

#Change the name of the fasta alignment file 
alignment=read.fasta("data/mothur_intermediate_files/Saanich.final.opti_mcc.unique_list.0.03.rep.fasta")
names(alignment)=rownames(taxonomy)

#remove OTUS not assigned to a domain
taxonomy=subset(taxonomy,!Domain=="unknown")

# Domain constrains
taxonomy[["Bacteria"]][taxonomy$Domain=="Bacteria"]=1
taxonomy[["Bacteria"]][taxonomy$Domain=="Archaea"]=0

# Phylum constrains
Phylum=unique(taxonomy$Phylum)
Phylum=subset(Phylum,!Phylum=="unknown_unclassified") #remove this factor

for (i in Phylum)
  {
  taxonomy[[as.character(i)]][taxonomy$Phylum==i]=1
  taxonomy[[as.character(i)]][!taxonomy$Phylum==i]=0
  }

Constrains=taxonomy[,c("Bacteria",as.character(Phylum))] #keep only the constrains 

#Convert to fasta file
sequences=list()
for (i in 1:dim(Constrains)[1]){sequences[[i]]=Constrains[i,]}
write.fasta(sequences, names=rownames(Constrains), file.out="My_outputs/Phylogenetic_Constrains.fasta", open = "w", nbchar = 60, as.string = FALSE)

#Prune the alignment to sequences with assigned domains
alignment=alignment[rownames(Constrains)]
write.fasta(alignment, names=names(alignment), file.out="My_outputs/Saanich.final.opti_mcc.unique_list.0.03.rep_Names_Modified.fasta", open = "w", nbchar = 60, as.string = FALSE)

# Part 1.2. Tree Reconstruction with FastTree
#################################################

# This part is NOT done in R, but directly on the terminal usifn the FastTree program 

# Run Fast Tree with taxonomic constrains
#Desktop/Programmes_Unix/FastTree -gtr -cat 20 -constraints Desktop/Recherche/En_cours/workshopMicrobiome/Contenu/MyFiles/Phylogenetic_Constrains.fasta -nt Desktop/Recherche/En_cours/workshopMicrobiome/Contenu/MyFiles/Saanich.final.opti_mcc.unique_list.0.03.rep_Names_Modified.fasta > Desktop/Recherche/En_cours/workshopMicrobiome/Contenu/MyFiles/Saanish_FastTree

# Run Fast Tree without constrains
#Desktop/Programmes_Unix/FastTree -gtr -cat 20 -nt Desktop/Recherche/En_cours/workshopMicrobiome/Contenu/MyFiles/Saanich.final.opti_mcc.unique_list.0.03.rep_Names_Modified.fasta > Desktop/Recherche/En_cours/workshopMicrobiome/Contenu/MyFiles/Saanish_FastTree_withoutConstrains

# Part 1.3. Tree Visualization 
#################################################

#Load the trees
Tree=read.tree('My_outputs/Saanish_FastTree')
TreeNoC=read.tree('My_outputs/Saanish_FastTree_withoutConstrains')

# Plot the tree with colours for Domains
Domains=unique(taxonomy$Domain)
paletteDomains=rainbow(length(Domains));names(paletteDomains)=Domains #define colours for taxonomic groups (domains here)
coloursDomains=paletteDomains[as.character(taxonomy$Domain)];names(coloursDomains)=rownames(taxonomy) #assign colors to each OTU depending on its taxonomic group 
pdf("My_outputs/Phylogenetic_tree_colouredby_Domains.pdf",width=15,height=15)
plot(Tree,type="fan",cex=.3,tip.color=coloursDomains[Tree$tip.label])
legend(1, 0, legend=names(paletteDomains),fill=paletteDomains, cex=2)
dev.off()

# Plot the tree with colours for Phylum
Phylums=unique(taxonomy$Phylum)
palettePhylums=rainbow(length(Phylums));names(palettePhylums)=Phylums #define colours for taxonomic groups (Phylums here)
coloursPhylums=palettePhylums[as.character(taxonomy$Phylum)];names(coloursPhylums)=rownames(taxonomy) #assign colors to each OTU depending on its taxonomic group 
pdf("My_outputs/Phylogenetic_tree_colouredby_Phylum.pdf",width=15,height=15)
plot(Tree,type="fan",cex=.3,tip.color=coloursPhylums[Tree$tip.label])
legend(1, 1, legend=names(palettePhylums),fill=palettePhylums, cex=1)
dev.off()

# check the impact of the topological contrains 
pdf("My_outputs/Phylogenetic_tree_colouredby_Phylum_NoConstrains.pdf",width=15,height=15)
plot(TreeNoC,type="fan",cex=.3,tip.color=coloursDomains[TreeNoC$tip.label])
legend(1, 1, legend=names(paletteDomains),fill=paletteDomains, cex=1)
dev.off()

# Part 1.4. Re-root the tree (Archaea Vs Bacteria)
#################################################

#First find the branch on which to root the tree 
Archaea=row.names(taxonomy)[taxonomy$Domain=="Archaea"]
MRCAnode=getMRCA(phy = Tree,tip = Archaea)

#re root
TreeNewRoot=root(phy=Tree, node=MRCAnode)

# Plot the tree with colours for Domains
pdf("My_outputs/Phylogenetic_tree_colouredby_Domains_NewRoot.pdf",width=15,height=15)
plot(TreeNewRoot,type="fan",cex=.3,tip.color=coloursDomains[Tree$tip.label])
legend(1, 0, legend=names(paletteDomains),fill=paletteDomains, cex=2)
dev.off()

################################################################################
# Part II: Diversity Analysis in R 
################################################################################



