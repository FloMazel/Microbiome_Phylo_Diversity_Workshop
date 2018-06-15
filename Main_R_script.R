################################################################################
# load and install packages
################################################################################

# if the named package is not installed then install it
library(seqinr)
library(ape)
library(vegan)
library(tidyverse)
library(phyloseq)
library(betapart)
library(abind)
library(tidyr)
library(Matrix)

#Put here your working directory (i.e you need to replace the location of the file to match the location of the downloaded folder on YOUR computer)
setwd("/Users/fmazel/Documents/GitHub/WorkingDirectory/")


################################################################################
# Part I: Construct the phylogenetic Tree
################################################################################

# Part 1.2. Construct the constrains file based on taxonomic information 
################################################################################

# load the taxonomic file
taxonomy.raw = read.table("data/Saanich_cruise72_mothur_OTU_taxonomy.taxonomy", sep="\t", header=TRUE, row.names=1)
taxonomy= taxonomy.raw %>% 
  select(-Size) %>% 
  separate(Taxonomy, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")

# load the alignment 
alignment=read.fasta("data/mothur_intermediate_files/Saanich.final.opti_mcc.unique_list.0.03.rep.fasta")
names(alignment)=rownames(taxonomy) #Change the name of the fasta alignment file 

#remove OTUS not assigned to a domain
taxonomy=subset(taxonomy,!Domain=="unknown")

# Domain constrains
taxonomy[["Bacteria"]][taxonomy$Domain=="Bacteria"]=1
taxonomy[["Bacteria"]][taxonomy$Domain=="Archaea"]=0

# Phylum constrains
Phylum=unique(taxonomy$Phylum)
Phylum=subset(Phylum,!Phylum%in%c("unknown_unclassified","Bacteria_unclassified", "Archaea_unclassified")) #remove this factor

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

# This part is NOT done in R, but directly on the terminal using the FastTree program 

# Run Fast Tree with taxonomic constrains
#Desktop/Programmes_Unix/FastTree -gtr -cat 20 -constraints /Users/fmazel/Documents/GitHub/WorkingDirectory/My_outputs/Phylogenetic_Constrains.fasta -nt /Users/fmazel/Documents/GitHub/WorkingDirectory/My_outputs/Saanich.final.opti_mcc.unique_list.0.03.rep_Names_Modified.fasta > /Users/fmazel/Documents/GitHub/WorkingDirectory/My_outputs/Saanish_FastTree

# Run Fast Tree without constrains
#Desktop/Programmes_Unix/FastTree -gtr -cat 20 -nt /Users/fmazel/Documents/GitHub/WorkingDirectory/My_outputs/Saanich.final.opti_mcc.unique_list.0.03.rep_Names_Modified.fasta > /Users/fmazel/Documents/GitHub/WorkingDirectory/My_outputs/Saanish_FastTree_withoutConstrains

# Part 1.3. Tree Visualization 
#################################################

# Load the trees
# Note: If you did not contruct the tree, replace 'My_outputs' by 'My_outputs_Back_Up'

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
plot(TreeNoC,type="fan",cex=.3,tip.color=coloursPhylums[TreeNoC$tip.label])
legend(1, 1, legend=names(palettePhylums),fill=palettePhylums, cex=1)
dev.off()

# Part 1.4. Re-root the tree (Archaea Vs Bacteria)
#################################################

#First find the branch on which to root the tree 
Archaea=row.names(taxonomy)[taxonomy$Domain=="Archaea"]
MRCAnode=getMRCA(phy = Tree,tip = Archaea)

#re root
TreeNewRoot=root(phy=Tree, node=MRCAnode,resolve.root = T)
#TreeNewRoot1=reroot(tree=Tree, node.number =MRCAnode)

TreeNewRoot
write.tree(TreeNewRoot,'My_outputs/Saanish_FastTreeRooted')

# Plot the tree with colours for Domains
pdf("My_outputs/Phylogenetic_tree_colouredby_Domains_NewRoot.pdf",width=15,height=15)
plot(TreeNewRoot,type="fan",cex=.3,tip.color=coloursDomains[Tree$tip.label])
legend(1, 0, legend=names(paletteDomains),fill=paletteDomains, cex=2)
dev.off()

################################################################################
# Part II: Diversity Analysis in R 
################################################################################

# 1. Import and clean up data

#Tree
# Note: If you did not contruct the tree, replace 'My_outputs' by 'My_outputs_Back_Up'

Tree=read.tree('My_outputs/Saanish_FastTreeRooted')

#OTU table
OTU = read.table("data/Saanich_cruise72_mothur_OTU_table.shared", sep="\t", header=TRUE, row.names=2)
OTU.clean = OTU %>%
  select(-label, -numOtus)

#taxonomy 
taxonomy = read.table("data/Saanich_cruise72_mothur_OTU_taxonomy.taxonomy", sep="\t", header=TRUE, row.names=1)
taxonomy.clean = taxonomy %>% 
  select(-Size) %>% 
  separate(Taxonomy, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")

#metadata
metadata = read.table("data/Saanich_cruise72_metadata.txt", sep="\t", header=TRUE, row.names=1)
metadata=metadata[,c("Depth_m","PO4_uM","SiO2_uM","NO3_uM","NH4_uM","CH4_nM" ,"Salinity_PSU")]
#colnames(metadata)
# Construct the phyloseq object 

#Define parts of the phyloseq object.
OTU.clean.physeq = otu_table(as.matrix(OTU.clean), taxa_are_rows=FALSE)
tax.clean.physeq = tax_table(as.matrix(taxonomy.clean))
metadata.physeq = sample_data(metadata)

phylogeny.physeq=phy_tree(Tree)

# apply(otu_table(mothur),1,sum) # Number of reads by sample
# At this point you can rarefy if needed

# 2. Taxonomic and phylogenetic  Beta-diversity

#Compute and compare metrics #NB -- Bray-Curtis and Jaccard are also implemeted in phyloseq; using vegdist FYI
BC=vegdist(otu_table(mothur),method = "bray")
Jaccard=vegdist(otu_table(mothur),method = "jac")
UniFracBeta=UniFrac(mothur)

#UniFracWBeta=UniFrac(mothur,weighted = T)
plot(Jaccard,UniFracBeta)

# Visualize Beta 
ordi = ordinate(mothur, "PCoA", "unifrac", weighted=F)
plot_ordination(mothur, ordi, color="Depth_m")
ggsave("My_outputs/PCoA_Unifrac.pdf",height = 4,width = 6)


ordi = ordinate(mothur, "PCoA", "bray", weighted=F)
plot_ordination(mothur, ordi, color="Depth_m")

# Test the relationship with environment using PERMANOVA
adonis(BC~Depth_m,data=data.frame(sample_data(mothur)))
adonis(Jaccard~Depth_m,data=data.frame(sample_data(mothur)))
adonis(UniFracBeta~Depth_m,data=data.frame(sample_data(mothur)))
adonis(UniFracBeta~Depth_m+NO3_uM,data=data.frame(sample_data(mothur)))

#the importance of phylogeny
adonis(Jaccard~NO3_uM,data=data.frame(sample_data(mothur)))
adonis(UniFracBeta~NO3_uM,data=data.frame(sample_data(mothur)))


# 3. BDTT
source("R functions/BDTT_functions.R")

#Recall the scale of divergence times
Hnodes=getHnodes(Tree) 
hist(Hnodes,n=150,xlim=c(0,.5))
#defines the slices
slices=c(seq(from=0,to=0.3,by=0.025)) 

#Run BDTT 
mat=t(as(otu_table(mothur), "matrix"))
MultipleBetaJac=BDTT(similarity_slices=slices,tree=Tree,sampleOTUs=mat,onlyBeta=T,metric="jac")
saveRDS(MultipleBetaJac,"My_outputs/Multiple_Resolution_Beta_Jaccard.RDS")  

#MultipleBetaBC=BDTT(similarity_slices=slices,tree=Tree,sampleOTUs=mat,onlyBeta=T,metric="bc")
#saveRDS(MultipleBetaBC,"My_outputs/Multiple_Resolution_Beta_BrayCurtis.RDS")  
#Because it is slow to run, we will not run bray curtis but directly load it from the backUp
MultipleBetaBC=readRDS("My_outputs/Multiple_Resolution_Beta_BrayCurtis.RDS")

# 3. Statistical tests: m predictors * n number of slices = a lot of models!
predictors=names(sample_data(mothur))
StatsRes=expand.grid(similarity_slices=as.character(slices),predictors=predictors,metric=c("Jac","BC"))
StatsRes[["F.Model"]]=StatsRes[["R2"]]=StatsRes[["Pr(>F)"]]=NA

for (i in as.character(slices))
{
  for (j in predictors) 
  {
   res=unlist(adonis(formula = MultipleBetaJac[i,,]~data.frame(sample_data(mothur))[,j])$aov.tab[1,c(6,5,4)])
   StatsRes[(StatsRes$metric=="Jac")&(StatsRes$predictors==j)&(StatsRes$similarity_slices==i),4:6]=res
   res=unlist(adonis(formula = MultipleBetaBC[i,,]~data.frame(sample_data(mothur))[,j])$aov.tab[1,c(6,5,4)])
   StatsRes[(StatsRes$metric=="BC")&(StatsRes$predictors==j)&(StatsRes$similarity_slices==i),4:6]=res
   }
}

# Plotting
ggplot(aes(y=F.Model,x=similarity_slices,colour=predictors,group=predictors),data=StatsRes)+geom_point()+geom_line()+facet_wrap(~metric)

ggplot(aes(y=R2,x=similarity_slices,colour=predictors,group=predictors),data=StatsRes)+geom_point()+geom_line()+facet_wrap(~metric)
ggsave("My_outputs/BDTT_Jaccard_BC.pdf",height = 7,width = 10)

#Put only the points that are significants
StatsRes$R2[StatsRes$`Pr(>F)`>.05]=NA
ggplot(aes(y=R2,x=similarity_slices,colour=predictors,group=predictors),data=StatsRes)+geom_point()+geom_line()+facet_wrap(~metric)

# Look more carefully at the effect of depth along the phylogenetic time scale
plot(as.matrix(dist(sample_data(mothur)[,"Depth_m"])),MultipleBetaBC["0",,])
plot(as.matrix(dist(sample_data(mothur)[,"Depth_m"])),MultipleBetaBC["0.125",,])


