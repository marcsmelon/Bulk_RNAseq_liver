.libPaths("/team_folders/cox_lab/R_libraries/4.0.2_Marcos")
setwd("/team_folders/cox_lab/Cox_Lab_Projects/BulkRNA_seq_cryoinjury")
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)


# SET THE DESIRED ORGANISM HERE
organism = "org.Dm.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

################
###SHAM vs 1 dpci
################

sham_1dpci_ensembl <- merge(zebrafish_annotation_geneinfo, sham_1dpci, all=FALSE)
sham_1dpci_ensembl_2 <- merge(sham_1dpci_ensembl, zb_human, all=FALSE)
write.table(sham_1dpci_ensembl_2, file="Sham_1dpci_ensembl_2.txt", col.names=TRUE,row.names=FALSE, quote=FALSE, sep="\t")
head(sham_1dpci_ensembl_1)
dim(sham_1dpci_ensembl_1)
dim(sham_1dpci_ensembl)
dim(sham_1dpci_ensembl)
head(sham_1dpci_ensembl)
View(sham_1dpci_ensembl)

# reading in data from deseq2
df = sham_1dpci
# we want the log2 fold change 
original_gene_list <- df$log2FoldChange
# name the vector
names(original_gene_list) <- df$LLgeneAbbrev
# omit any NA values 
gene_list<-na.omit(original_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)


gse1<- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

require(DOSE)
dotplot(gse1, showCategory=30, split=".sign") + facet_grid(.~.sign)
ridgeplot(gse) + labs(x = "enrichment distribution")
gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 10)




# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = df[df$LLgeneAbbrev %in% dedup_ids$SYMBOL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

kegg_organism = "dme"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")
dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)

################
###SHAM vs 3 dpci
################

# reading in data from deseq2
df1 = read.csv("sham_3dpci.csv", header=TRUE)
# we want the log2 fold change 
original_gene_list1 <- df1$log2FoldChange
# name the vector
names(original_gene_list1) <- df1$geneName
# omit any NA values 
gene_list1<-na.omit(original_gene_list1)
# sort the list in decreasing order (required for clusterProfiler)
gene_list1 = sort(gene_list1, decreasing = TRUE)


gse2 <- gseGO(geneList=gene_list1, 
             ont ="ALL", 
             keyType = "GENENAME", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

require(DOSE)
dotplot(gse2, showCategory=20, split=".sign") + facet_grid(.~.sign)
ridgeplot(gse) + labs(x = "enrichment distribution")
gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 10)



################
###SHAM vs 7 dpci
################

# reading in data from deseq2
df2 = read.csv("sham_7dpci.csv", header=TRUE)
# we want the log2 fold change 
original_gene_list2 <- df2$log2FoldChange
# name the vector
names(original_gene_list2) <- df2$geneName
# omit any NA values 
gene_list2<-na.omit(original_gene_list2)
# sort the list in decreasing order (required for clusterProfiler)
gene_list2 = sort(gene_list2, decreasing = TRUE)


gse3 <- gseGO(geneList=gene_list2, 
             ont ="ALL", 
             keyType = "GENENAME", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

require(DOSE)
dotplot(gse3, showCategory=20, split=".sign") + facet_grid(.~.sign)
ridgeplot(gse) + labs(x = "enrichment distribution")
gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 10)


# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(original_gene_list2), fromType = "G", toType = "ENTREZID", OrgDb=organism)
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

